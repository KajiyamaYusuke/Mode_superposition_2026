#!/usr/bin/env python3
"""Read-only area.dat analysis; never imports or invokes simulation drivers."""

import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[2]
CASES = {
    "R0": ROOT / "output/kanaya_ablation/R0_modern",
    "S2": ROOT / "output/kanaya_ablation/S2_legacy_vt0175",
    "S3": ROOT / "output/kanaya_ablation_next/S3_legacy_upstream_circuit",
}
DEFAULT_OUTPUT = ROOT / "output/kanaya_ablation_next/area_waveform_comparison_010_020"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimum_area(data):
    """NaN marks an undefined section, not zero area or a closed section."""
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError("expected step followed by section areas")
    areas = data[:, 1:]
    if np.isinf(areas).any() or (areas < 0).any():
        raise ValueError("infinite or negative section area")
    if not np.isfinite(areas).any(axis=1).all():
        raise ValueError("no valid section in a sampled row")
    result = np.nanmin(areas, axis=1)
    # Independent row-wise reduction, including the invalid-section handling.
    reference = np.array([min(row[np.isfinite(row)]) for row in areas])
    np.testing.assert_array_equal(result, reference)
    return result


def spectrum(area, sample_dt):
    window = np.hanning(len(area))
    amplitudes = 2 * np.abs(np.fft.rfft((area - area.mean()) * window)) / window.sum()
    amplitudes[0] = 0
    amplitudes[-1] *= 0.5 if len(area) % 2 == 0 else 1
    frequency = np.fft.rfftfreq(len(area), sample_dt)
    peak = int(np.argmax(amplitudes[1:]) + 1)
    return frequency, amplitudes, peak


def opening_peaks(area, time, prominence_fraction):
    peaks, properties = find_peaks(area, prominence=prominence_fraction * np.ptp(area))
    if len(peaks) < 2:
        raise ValueError("fewer than two prominent opening peaks")
    periods = np.diff(time[peaks])
    frequency = 1 / periods.mean()
    np.testing.assert_allclose(frequency, (len(peaks) - 1) / (time[peaks[-1]] - time[peaks[0]]))
    return peaks, properties, periods, float(frequency)


def analyze(label, directory, start, end, prominence_fraction):
    area_path = directory / "output/dat/area.dat"
    manifest_path = directory / "manifest.json"
    sources = [area_path, manifest_path]
    before = {str(p): {"sha256": sha(p), "bytes": p.stat().st_size,
                       "mtime_ns": p.stat().st_mtime_ns} for p in sources}
    manifest = json.loads(manifest_path.read_text())
    dt = float(manifest["effective_parameters"]["dt_s"])
    raw = np.loadtxt(area_path, ndmin=2)
    steps = raw[:, 0]
    if not np.isfinite(steps).all() or not np.equal(steps, np.round(steps)).all():
        raise ValueError("step column must contain finite integers")
    if not (np.diff(steps) > 0).all():
        raise ValueError("step numbers must be strictly increasing")
    # Compare in step units so decimal time roundoff cannot drop an endpoint.
    mask = (steps >= start / dt - 1e-8) & (steps <= end / dt + 1e-8)
    data = raw[mask]
    if len(data) < 3:
        raise ValueError("not enough samples in the requested interval")
    stride = np.diff(data[:, 0])
    if not np.equal(stride, stride[0]).all():
        raise ValueError("FFT requires uniformly sampled data")
    time = data[:, 0] * dt
    area = minimum_area(data)
    if np.ptp(area) == 0:
        raise ValueError("constant area has no identifiable vibration frequency")
    sample_dt = float(stride[0] * dt)
    frequency, amplitudes, fft_index = spectrum(area, sample_dt)
    peaks, properties, periods, peak_frequency = opening_peaks(area, time, prominence_fraction)
    rectangular = abs(np.fft.rfft(area - area.mean()))
    rectangular_index = int(np.argmax(rectangular[1:]) + 1)
    sensitivity = []
    for fraction in [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]:
        alternative, _, _, f = opening_peaks(area, time, fraction)
        sensitivity.append(dict(case=label, prominence_fraction=fraction,
                                peaks=len(alternative), mean_frequency_Hz=f,
                                same_peaks_as_selected=bool(np.array_equal(alternative, peaks))))
    mean = math.fsum(area) / len(area)
    np.testing.assert_allclose(mean, area.mean(), rtol=1e-14, atol=1e-14)
    summary = dict(case=label, maximum_mm2=float(area.max()), minimum_mm2=float(area.min()),
                   mean_mm2=float(mean), FFT_peak_Hz=float(frequency[fft_index]),
                   opening_peak_mean_frequency_Hz=peak_frequency,
                   mean_period_ms=float(periods.mean() * 1000),
                   minimum_period_ms=float(periods.min() * 1000),
                   maximum_period_ms=float(periods.max() * 1000),
                   peak_count=len(peaks), interval_count=len(periods), samples=len(area),
                   actual_start_s=float(time[0]), actual_end_s=float(time[-1]),
                   sample_dt_s=sample_dt, FFT_bin_spacing_Hz=float(frequency[1]),
                   maximum_time_s=float(time[np.argmax(area)]), minimum_time_s=float(time[np.argmin(area)]),
                   invalid_section_entries=int(np.isnan(data[:, 1:]).sum()),
                   rectangular_FFT_peak_Hz=float(frequency[rectangular_index]),
                   prominence_fraction=prominence_fraction)
    waveform = pd.DataFrame(dict(case=label, step=data[:, 0].astype(int), time_s=time,
                                 Amin_mm2=area, invalid_sections=np.isnan(data[:, 1:]).sum(axis=1)))
    peak_table = pd.DataFrame(dict(case=label, peak_number=np.arange(1, len(peaks) + 1),
                                  step=data[peaks, 0].astype(int), time_s=time[peaks],
                                  Amin_mm2=area[peaks], prominence_mm2=properties["prominences"]))
    interval_table = pd.DataFrame(dict(case=label, interval_number=np.arange(1, len(periods) + 1),
                                      start_s=time[peaks[:-1]], end_s=time[peaks[1:]],
                                      period_s=periods, instantaneous_frequency_Hz=1 / periods))
    fft_table = pd.DataFrame(dict(case=label, frequency_Hz=frequency, amplitude_mm2=amplitudes))
    return summary, waveform, peak_table, interval_table, fft_table, sensitivity, before


def plots(output, summary, waves, peaks, spectra):
    os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="area-comparison-mpl-"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"R0": "#0072B2", "S2": "#D55E00", "S3": "#009E73"}
    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True, sharey=True)
    for ax, row in zip(axes, summary.itertuples()):
        wave = waves[waves.case == row.case]
        peak = peaks[peaks.case == row.case]
        ax.plot(wave.time_s, wave.Amin_mm2, color=colors[row.case], lw=1, label="Amin")
        ax.plot(peak.time_s, peak.Amin_mm2, "o", color="black", ms=3, label="opening peaks")
        ax.axhline(row.mean_mm2, color=colors[row.case], ls="--", alpha=.6, label="mean")
        ax.set_title(f"{row.case}: max {row.maximum_mm2:.4f}, min {row.minimum_mm2:.5f}, "
                     f"mean {row.mean_mm2:.4f} mm$^2$; FFT {row.FFT_peak_Hz:.0f} Hz, "
                     f"peak intervals {row.opening_peak_mean_frequency_Hz:.2f} Hz", fontsize=10)
        ax.set_ylabel("Minimum area [mm$^2$]")
        ax.grid(alpha=.25)
        ax.legend(loc="upper left", ncol=3, fontsize=8)
    axes[-1].set_xlabel("Time [s]")
    axes[-1].set_xlim(.10, .20)
    axes[-1].set_ylim(bottom=0)
    fig.tight_layout()
    for suffix in ["png", "svg"]:
        fig.savefig(output / f"area_waveforms.{suffix}", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    for row in summary.itertuples():
        spec = spectra[(spectra.case == row.case) & (spectra.frequency_Hz <= 1000)]
        ax.plot(spec.frequency_Hz, spec.amplitude_mm2, color=colors[row.case],
                label=f"{row.case}: {row.FFT_peak_Hz:.0f} Hz", marker=".", ms=3)
    ax.set(xlabel="Frequency [Hz]", ylabel="Hann-window amplitude [mm$^2$]",
           title="Minimum-area spectrum (mean removed; no zero padding)", xlim=(0, 1000))
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "area_fft.png", dpi=180)
    plt.close(fig)


def report(output, summary, cases):
    lines = ["# R0・S2・S3 最狭断面積波形の比較", "",
             "既存の `area.dat` のみを波形入力として後処理した。連成計算・形状再構成・新規物理計算は実行していない。", "",
             "## 比較結果", "",
             "| ケース | 最大値 mm² | 最小値 mm² | 平均値 mm² | FFTピーク Hz | 開口ピーク間隔による平均周波数 Hz | ピーク数／間隔数 |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for r in summary.itertuples():
        lines.append(f"| {r.case} | {r.maximum_mm2:.5f} | {r.minimum_mm2:.7f} | {r.mean_mm2:.6f} | "
                     f"{r.FFT_peak_Hz:.2f} | {r.opening_peak_mean_frequency_Hz:.4f} | {r.peak_count}／{r.interval_count} |")
    lines += ["", "![最狭断面積波形と主開口ピーク](area_waveforms.png)", "",
              "![面積波形のFFT](area_fft.png)", "", "## 定義・処理条件", "",
              "- 各行の1列目はstep、以降は各断面の面積（mm²）。同じ時刻の有限な断面値の最小を `Amin(t)` とした。最狭断面の位置は時刻ごとに変わり得る。R0のNaNは幾何的に無効な断面であり、0への置換や時間方向の補間はしない。無効断面を除いても全時刻に有効値があることを確認した。",
              "- 元manifestのdt=10 µs、area.datの出力間隔は5step=50 µs（20 kHz）。指定区間0.10〜0.20 sに実在するデータは各ケースとも **0.10000〜0.19995 sの2,000点**。0.20000 sは保存されていないため、外挿・補完しない。平均はこの2,000点の算術平均。",
              "- FFTはAminの平均を除去し、対称Hann窓を適用した一側FFTの非DC最大ビン。ゼロ詰め・ピーク補間なし。ビン間隔は `1/(N Δt)=10 Hz`。窓なしFFTでも最大ビンは3ケースとも同じだった。FFT値の整数差を10 Hzより細かい周波数差の証拠にしない。",
              "- 開口ピークは未平滑化Aminから `scipy.signal.find_peaks` で抽出。prominenceは各ケースの `(最大−最小)×0.20`。FFT由来の最小距離や周期を設定せず、独立に抽出した。区間端の未完ピークは含めない。平均周波数は **`1 / mean(t_peak[i+1]−t_peak[i])`**。瞬時周波数の算術平均とは区別する。",
              "- prominenceを20/25/30/40/50%へ変えても、3ケースとも同一ピーク・同一周波数になった。R0では10/15%に下げると副ピークまで計数されるため、全局所極大を一周期とみなす周波数は採用しない。しきい値感度は `peak_sensitivity.csv` に保存した。",
              "- 過去レポートの高精度・毎step診断とはサンプリングと保存精度が異なる。今回は指定どおり丸め済みarea.datだけから統計量を求めたため、過去の最小値・平均値と末尾が異なり得る。",
              "- これらは面積波形の振動周波数であり、閉鎖を必要とするCQや閉鎖由来f0ではない。", "",
              "## ピーク間隔の確認", "",
              "| ケース | 平均周期 ms | 最短〜最長周期 ms |",
              "|---|---:|---:|"]
    for r in summary.itertuples():
        lines.append(f"| {r.case} | {r.mean_period_ms:.6f} | {r.minimum_period_ms:.5f}〜{r.maximum_period_ms:.5f} |")
    lines += ["", "R0→S2→S3で主振動周波数は上昇し、開口ピーク間隔の値はそれぞれ対応するFFT最大ビンの近傍にある。S3はS2より面積最大値・最小値・平均値がいずれも大きい。3ケースとも、この保存波形の最小値は公式閉鎖しきい値0.01 mm²より大きい。", "",
              "## 入力・成果物", ""]
    for label, directory in cases.items():
        lines.append(f"- {label}: `{directory / 'output/dat/area.dat'}`")
    lines += ["", "`comparison.csv`：比較表・極値時刻・サンプル条件。`waveforms.csv`：全6,000点のAmin。`opening_peaks.csv`：採用ピーク時刻・面積・prominence。`peak_intervals.csv`：全隣接周期と瞬時周波数。`fft_spectra.csv`：全FFTビン。`provenance.json`：入力・解析スクリプトSHA256と入力不変検証。", "",
              "再実行：別の未使用出力先を指定して `python3 -B tools/kanaya_ablation/compare_area_waveforms.py --output <new-directory>`。既存出力先への上書きは拒否する。", ""]
    (output / "report.md").write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--r0", type=Path, default=CASES["R0"], help="R0 case directory")
    parser.add_argument("--s2", type=Path, default=CASES["S2"], help="S2 case directory")
    parser.add_argument("--s3", type=Path, default=CASES["S3"], help="S3 case directory")
    args = parser.parse_args()
    cases = {"R0": args.r0.resolve(), "S2": args.s2.resolve(), "S3": args.s3.resolve()}
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    results = [analyze(label, directory, .10, .20, .20) for label, directory in cases.items()]
    summary = pd.DataFrame([r[0] for r in results])
    waves, peaks, intervals, spectra = [pd.concat([r[i] for r in results], ignore_index=True) for i in range(1, 5)]
    sensitivity = pd.DataFrame([row for r in results for row in r[5]])
    assert sensitivity[sensitivity.prominence_fraction >= .20].same_peaks_as_selected.all()
    assert np.array_equal(summary.FFT_peak_Hz, summary.rectangular_FFT_peak_Hz)
    output.mkdir(parents=True, exist_ok=False)
    for name, frame in [("comparison", summary), ("waveforms", waves), ("opening_peaks", peaks),
                        ("peak_intervals", intervals), ("fft_spectra", spectra), ("peak_sensitivity", sensitivity)]:
        frame.to_csv(output / f"{name}.csv", index=False)
    plots(output, summary, waves, peaks, spectra)
    report(output, summary, cases)
    inputs = {p: info for r in results for p, info in r[6].items()}
    for path, info in inputs.items():
        p = Path(path)
        assert info == {"sha256": sha(p), "bytes": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
    provenance = dict(passed=True, input_files_unchanged=True, coupled_simulation_executed=False,
                      inputs=inputs, script=str(Path(__file__).resolve()), script_sha256=sha(Path(__file__)),
                      requested_window_s=[.10, .20], peak_prominence_fraction=.20,
                      FFT="demean, symmetric Hann, non-DC max bin, no padding/interpolation",
                      peak_frequency="reciprocal of mean adjacent prominent opening peak interval",
                      all_nan_rows=0, infinite_areas=0, negative_areas=0,
                      peak_sensitivity_20_to_50_percent_identical=True,
                      rectangular_and_Hann_FFT_peak_bins_identical=True,
                      numpy_version=np.__version__)
    (output / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n")
    print(summary[["case", "maximum_mm2", "minimum_mm2", "mean_mm2", "FFT_peak_Hz",
                   "opening_peak_mean_frequency_Hz", "peak_count"]].to_string(index=False))
    print("Results:", output)


if __name__ == "__main__":
    main()

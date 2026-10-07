import subprocess
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import scienceplots
from steinecke_classification import classify_locking, one_sided_spectrum

# 論文用のスタイル設定
plt.style.use(['science', 'ieee', 'no-latex'])

# =========================
# ユーザー設定環境
# =========================
EXECUTABLE = "/home/kajiyama/code/simulation-main/build/simulation"
PARAM_FILE = "/home/kajiyama/code/simulation-main/input/param.txt"
OUTPUT_DATA = "/home/kajiyama/code/simulation-main/output/dat/airflow_vt.dat"
DISP_DATA = "/home/kajiyama/code/simulation-main/output/dat/displace.dat"
AREA_DATA = "/home/kajiyama/code/simulation-main/output/dat/area.dat"
OUTPUT_DIR = Path("/home/kajiyama/code/simulation-main/output/png")
CSV_DIR = Path("/home/kajiyama/code/simulation-main/output/csv")
CLASSIFICATION_PNG_DIR = OUTPUT_DIR / "classification"
CLASSIFICATION_CSV_DIR = CSV_DIR / "classification"

# スイープする圧力のリスト (Pa)
pressure_list = [600, 1000, 1400]
# 左右に共通して設定するモード減衰率
zeta_list = [0.1]

# 解析設定
sim_dt = 1.0e-5
output_interval = 5  
dt = sim_dt * output_interval
fs = 1.0 / dt         
t_start = 0.0 
t_end   = 0.5        
flow_plot_start = 0.2
flow_plot_end = 0.3

# Steinecke--Herzel classification helper thresholds. These are numerical
# helper rules, not thresholds specified by the original paper.
classification_duration = 0.15
classification_prominence_fraction = 0.05
classification_max_period = 12
classification_min_peaks = 8
classification_amplitude_rtol = 0.02
classification_interval_rtol = 0.02
classification_count_agreement = 0.8
classification_subharmonic_ratio = 0.05

M3_PER_S_TO_L_PER_S = 1000.0


def flow_m3s_to_ls(flow_m3s):
    """Convert volumetric flow from m^3/s (solver output) to L/s."""
    return np.asarray(flow_m3s) * M3_PER_S_TO_L_PER_S


def minimum_glottal_area(area_data):
    """Return the minimum section area [mm^2] for every saved step."""
    area_data = np.atleast_2d(area_data)
    if area_data.shape[1] < 2:
        raise ValueError("area.dat must contain a step and at least one area")
    return np.min(area_data[:, 1:], axis=1)


def load_complete_rows(path, minimum_columns):
    """Load numeric rows while ignoring an incomplete trailing output row."""
    rows = []
    expected_columns = None
    with open(path, 'r', encoding='utf-8') as stream:
        for line in stream:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            try:
                row = [float(value) for value in stripped.split()]
            except ValueError:
                continue
            if len(row) < minimum_columns:
                continue
            if expected_columns is None:
                expected_columns = len(row)
            if len(row) == expected_columns:
                rows.append(row)
    if not rows:
        raise ValueError(f"{path} contains no complete numeric rows")
    return np.asarray(rows)


def condition_tag(pressure_val, zeta_val):
    zeta_text = f"{zeta_val:g}".replace('.', 'p')
    return f"Ps_{pressure_val}Pa_zeta_{zeta_text}"


def _write_rows(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def classify_and_save(pressure_val, zeta_val):
    """Classify a steady response and save every diagnostic used by the label."""
    data_disp = load_complete_rows(DISP_DATA, 3)
    time = data_disp[:, 0]
    left = -data_disp[:, 1]
    right = data_disp[:, 2]
    analysis_end = min(t_end, float(time[-1]))
    analysis_start = max(t_start, analysis_end - classification_duration)
    result = classify_locking(
        time, left, right,
        transient_end=analysis_start,
        prominence_fraction=classification_prominence_fraction,
        max_period=classification_max_period,
        min_peaks=classification_min_peaks,
        amplitude_rtol=classification_amplitude_rtol,
        interval_rtol=classification_interval_rtol,
        count_agreement=classification_count_agreement,
    )
    if analysis_end - analysis_start < 0.95 * classification_duration:
        result['label'] = 'insufficient_data'
        result.setdefault('consistency_issues', []).append('analysis_window_too_short')
    tag = condition_tag(pressure_val, zeta_val)
    CLASSIFICATION_CSV_DIR.mkdir(parents=True, exist_ok=True)
    CLASSIFICATION_PNG_DIR.mkdir(parents=True, exist_ok=True)

    peak_rows = []
    for side, peaks in (
        ('right', result.get('peaks_right')),
        ('left', result.get('peaks_left')),
    ):
        if peaks is None:
            continue
        peak_rows.extend(
            {'side': side, 'peak_number': number + 1, 'sample_index': int(index),
             'time_s': peak_time, 'displacement_mm': value}
            for number, (index, peak_time, value) in enumerate(
                zip(peaks.indices, peaks.times, peaks.values)
            )
        )
    _write_rows(
        CLASSIFICATION_CSV_DIR / f"peaks_{tag}.csv",
        ['side', 'peak_number', 'sample_index', 'time_s', 'displacement_mm'],
        peak_rows,
    )

    windows = np.asarray(result.get('count_windows', np.empty((0, 4))))
    _write_rows(
        CLASSIFICATION_CSV_DIR / f"total_cycle_windows_{tag}.csv",
        ['window_start_s', 'window_end_s', 'n_right', 'n_left'],
        [dict(zip(['window_start_s', 'window_end_s', 'n_right', 'n_left'], row))
         for row in windows],
    )
    error_rows = []
    for side, errors in (
        ('right', result.get('period_errors_right', {})),
        ('left', result.get('period_errors_left', {})),
    ):
        for period, (amplitude_error, interval_error) in errors.items():
            error_rows.append({
                'side': side, 'candidate_period': period,
                'normalized_amplitude_rmse': amplitude_error,
                'normalized_interval_rmse': interval_error,
            })
    _write_rows(
        CLASSIFICATION_CSV_DIR / f"period_errors_{tag}.csv",
        ['side', 'candidate_period', 'normalized_amplitude_rmse',
         'normalized_interval_rmse'], error_rows,
    )

    analysis_time = result.get('analysis_time', np.array([]))
    analysis_left = result.get('left_signal', np.array([]))
    analysis_right = result.get('right_signal', np.array([]))
    if len(analysis_time) >= 3:
        frequency, spectrum_left = one_sided_spectrum(analysis_time, analysis_left)
        _, spectrum_right = one_sided_spectrum(analysis_time, analysis_right)
    else:
        frequency = spectrum_left = spectrum_right = np.array([])
    subharmonic_ratio = np.nan
    subharmonic_clear = ''
    peaks_right_for_fft = result.get('peaks_right')
    if peaks_right_for_fft is not None and len(peaks_right_for_fft.times) >= 3 and len(frequency):
        base_frequency = 1.0 / float(np.median(np.diff(peaks_right_for_fft.times)))
        base_index = int(np.argmin(np.abs(frequency - base_frequency)))
        half_index = int(np.argmin(np.abs(frequency - 0.5 * base_frequency)))
        subharmonic_ratio = float(
            spectrum_right[half_index]
            / max(spectrum_right[base_index], np.finfo(float).eps)
        )
        subharmonic_clear = subharmonic_ratio >= classification_subharmonic_ratio
        if result.get('period_right') == 2 and not subharmonic_clear:
            issues = result.setdefault('consistency_issues', [])
            issues.append('subharmonic_not_clear')
            result['label'] = 'uncertain'
    _write_rows(
        CLASSIFICATION_CSV_DIR / f"spectrum_{tag}.csv",
        ['frequency_hz', 'left_amplitude', 'right_amplitude'],
        [{'frequency_hz': f, 'left_amplitude': sl, 'right_amplitude': sr}
         for f, sl, sr in zip(frequency, spectrum_left, spectrum_right)],
    )

    map_x = result.get('next_maximum_x', np.array([]))
    map_y = result.get('next_maximum_y', np.array([]))
    _write_rows(
        CLASSIFICATION_CSV_DIR / f"next_maximum_map_{tag}.csv",
        ['x_i_max_mm', 'x_next_max_mm'],
        [{'x_i_max_mm': x, 'x_next_max_mm': y} for x, y in zip(map_x, map_y)],
    )

    fig, axes = plt.subplots(2, 2, figsize=(9, 7), dpi=150)
    axes[0, 0].plot(analysis_time, analysis_left, label='Left', linewidth=0.8)
    axes[0, 0].plot(analysis_time, analysis_right, label='Right', linewidth=0.8)
    for side, peaks, color in (
        ('Left', result.get('peaks_left'), '#0072B2'),
        ('Right', result.get('peaks_right'), '#D55E00'),
    ):
        if peaks is not None:
            axes[0, 0].scatter(peaks.times, peaks.values, s=12, color=color, label=f'{side} peaks')
    axes[0, 0].set(xlabel='Time [s]', ylabel='Displacement [mm]', title='Steady waveform and peaks')
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].plot(analysis_left, analysis_right, linewidth=0.7)
    axes[0, 1].set(xlabel='Left displacement [mm]', ylabel='Right displacement [mm]', title='Left-right phase portrait')
    axes[1, 0].scatter(map_x, map_y, s=15)
    axes[1, 0].set(xlabel=r'$X_i^{max}$ [mm]', ylabel=r'$X_{i+1}^{max}$ [mm]', title='Right next-maximum map')
    if len(frequency):
        scale_l = max(float(np.max(spectrum_left)), np.finfo(float).eps)
        scale_r = max(float(np.max(spectrum_right)), np.finfo(float).eps)
        axes[1, 1].semilogy(frequency, spectrum_left / scale_l, label='Left')
        axes[1, 1].semilogy(frequency, spectrum_right / scale_r, label='Right')
    axes[1, 1].set(xlabel='Frequency [Hz]', ylabel='Normalized amplitude', title='One-sided spectra', xlim=(0, 1000))
    axes[1, 1].legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(True, linestyle='--', alpha=0.35)
        axis.tick_params(direction='in')
    base_label = result.get('base_label', result['label'])
    fig.suptitle(
        f"{tag}: classification={result['label']} (count label={base_label})",
        fontsize=13,
    )
    fig.tight_layout()
    figure_path = CLASSIFICATION_PNG_DIR / f"classification_{tag}.png"
    fig.savefig(figure_path, dpi=300)
    plt.close(fig)

    peaks_left = result.get('peaks_left')
    peaks_right = result.get('peaks_right')
    resolution = float(frequency[1] - frequency[0]) if len(frequency) > 1 else np.nan
    record = {
        'Q': '',  # This full simulation currently has no Steinecke Q parameter.
        'pressure_pa': pressure_val,
        'zeta': zeta_val,
        'analysis_start_s': analysis_start,
        'analysis_end_s': analysis_end,
        'label': result['label'],
        'count_label': result.get('base_label', ''),
        'n_right': result.get('n_right', ''),
        'n_left': result.get('n_left', ''),
        'period_right': result.get('period_right', ''),
        'period_left': result.get('period_left', ''),
        'period_errors_right': json.dumps(result.get('period_errors_right', {}), sort_keys=True),
        'period_errors_left': json.dumps(result.get('period_errors_left', {}), sort_keys=True),
        'total_cycle_s': result.get('T_total', ''),
        'count_agreement': result.get('count_agreement', ''),
        'number_of_peaks_right': 0 if peaks_right is None else len(peaks_right.times),
        'number_of_peaks_left': 0 if peaks_left is None else len(peaks_left.times),
        'next_maximum_cluster_count': result.get('map_cluster_count', ''),
        'fft_resolution_hz': resolution,
        'subharmonic_f_over_2_to_f_ratio': subharmonic_ratio,
        'subharmonic_check_passed': subharmonic_clear,
        'consistency_issues': ';'.join(result.get('consistency_issues', [])),
        'classification_thresholds': json.dumps({
            **result['thresholds'],
            'analysis_duration_s': classification_duration,
            'subharmonic_ratio_required': classification_subharmonic_ratio,
        }, sort_keys=True),
        'initial_condition_identifier': 'simulation_default',
    }
    print(f"[分類] {tag}: {record['label']} (count={record['count_label']})")
    return record


def write_classification_summary(records):
    if not records:
        return
    _write_rows(CSV_DIR / 'steinecke_classification_summary.csv', list(records[0]), records)


def save_full_history(pressure_val, zeta_val):
    """Save whole-run flow, minimum-area, and displacement histories."""
    data_flow = load_complete_rows(OUTPUT_DATA, 2)
    data_area = load_complete_rows(AREA_DATA, 2)
    data_disp = load_complete_rows(DISP_DATA, 3)

    flow_time = data_flow[:, 0] * sim_dt
    area_time = data_area[:, 0] * sim_dt
    disp_time = data_disp[:, 0]
    flow_ls = flow_m3s_to_ls(data_flow[:, 1])
    min_area_mm2 = minimum_glottal_area(data_area)
    left_displacement = -data_disp[:, 1]
    right_displacement = data_disp[:, 2]

    fig, axes = plt.subplots(3, 1, figsize=(8, 8), dpi=150, sharex=True)
    axes[0].plot(flow_time, flow_ls, color='black', linewidth=0.8)
    axes[0].set_ylabel("Airflow [L/s]", fontsize=13)

    axes[1].plot(area_time, min_area_mm2, color='#009E73', linewidth=0.8)
    axes[1].set_ylabel("Minimum glottal\narea [mm$^2$]", fontsize=13)

    axes[2].plot(
        disp_time, left_displacement, label='Left VF', color='#0072B2',
        linewidth=0.8,
    )
    axes[2].plot(
        disp_time, right_displacement, label='Right VF', color='#D55E00',
        linewidth=0.8, linestyle='--',
    )
    axes[2].set_ylabel("Displacement [mm]", fontsize=13)
    axes[2].set_xlabel("Time [s]", fontsize=14)
    axes[2].legend(loc='upper right', fontsize=11)

    final_time = max(flow_time[-1], area_time[-1], disp_time[-1])
    for ax in axes:
        ax.set_xlim(0.0, final_time)
        ax.tick_params(direction='in')
        ax.grid(True, linestyle='--', alpha=0.4)
    axes[0].set_title(
        f"Whole-run response (Ps = {pressure_val} Pa, zeta = {zeta_val:g})",
        fontsize=15,
    )
    fig.tight_layout()

    save_filename = OUTPUT_DIR / f"full_history_{condition_tag(pressure_val, zeta_val)}.png"
    fig.savefig(save_filename, dpi=300)
    plt.close(fig)
    print(f"[出力完了] 全時間履歴を保存しました: {save_filename}")

# =========================
# 関数: paramファイルの圧力書き換え
# =========================
def update_param_file(filepath, new_pressure, new_zeta):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    replacements = {
        "# damping coefficient zetaL": new_zeta,
        "# damping coefficient zetaR": new_zeta,
        "# subglottal pressure Ps (Pa)": new_pressure,
    }
    replaced = set()
    for i, line in enumerate(lines):
        for marker, value in replacements.items():
            if marker in line:
                lines[i + 1] = f"{value:g}\n"
                replaced.add(marker)

    missing = set(replacements) - replaced
    if missing:
        raise ValueError(f"param file is missing markers: {sorted(missing)}")
            
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(lines)
    print(
        f"[設定更新] 声門下圧={new_pressure} Pa, "
        f"zetaL=zetaR={new_zeta:g}"
    )

def save_flow_waveform(pressure_val, zeta_val):
    data_flow = load_complete_rows(OUTPUT_DATA, 2)
    time = data_flow[:, 0] * sim_dt
    flow = flow_m3s_to_ls(data_flow[:, 1])

    mask = (time >= flow_plot_start) & (time <= flow_plot_end)
    if not np.any(mask):
        print(
            f"[スキップ] 流量波形 {flow_plot_start:.1f}-{flow_plot_end:.1f} s のデータがありません。"
            f" 現在の最終時刻は {time[-1]:.5f} s です。"
        )
        return

    fig, ax = plt.subplots(figsize=(7, 3), dpi=150)
    ax.plot(time[mask], flow[mask]*1e3, color='black', linewidth=1.0)
    ax.set_xlabel("Time [s]", fontsize=14)
    ax.set_ylabel("Airflow [mL/s]", fontsize=14)
    ax.set_title(
        f"Airflow waveform (Ps = {pressure_val} Pa, zeta = {zeta_val:g})",
        fontsize=15,
    )
    ax.tick_params(direction='in')
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    save_filename = OUTPUT_DIR / (
        f"airflow_waveform_{flow_plot_start:.1f}_{flow_plot_end:.1f}s_"
        f"{condition_tag(pressure_val, zeta_val)}.png"
    )
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 流量波形を保存しました: {save_filename}")

# =========================
# 関数: シミュレーション結果の可視化
# =========================
def analyze_and_plot(pressure_val, zeta_val):
    # 2つのデータを読み込む
    data_flow = load_complete_rows(OUTPUT_DATA, 2)
    data_disp = load_complete_rows(DISP_DATA, 3)
    
    # 各ファイルに保存された時間・step列を使用する。
    time = data_flow[:, 0] * sim_dt
    flow_data = flow_m3s_to_ls(data_flow[:, 1])
    
    # displace.datから変位を取得 (1列目=L, 2列目=R)
    # 左声帯は符号が右声帯と逆なので、重ねて比較しやすいように反転する
    disp_time = data_disp[:, 0]
    x1l_data = -data_disp[:, 1]
    x1r_data = data_disp[:, 2]

    # 指定区間の切り出し
    flow_mask = (time >= t_start) & (time <= t_end)
    disp_mask = (disp_time >= t_start) & (disp_time <= t_end)
    valid_time = time[flow_mask]
    valid_flow = flow_data[flow_mask]
    valid_disp_time = disp_time[disp_mask]
    valid_x1l = x1l_data[disp_mask]
    valid_x1r = x1r_data[disp_mask]
    
    valid_flow_ac = valid_flow - np.mean(valid_flow)

    # スペクトログラム計算
    nperseg = 1024 
    noverlap = nperseg * 3 // 4 
    f, t_spec, Sxx = signal.spectrogram(valid_flow_ac, fs=fs, window='hann',
                                        nperseg=nperseg, noverlap=noverlap)
    t_spec = t_spec + t_start

    Sxx_linear = np.sqrt(Sxx)
    max_Sxx = np.max(Sxx_linear) if np.max(Sxx_linear) > 0 else 1.0
    Sxx_db = 20 * np.log10(Sxx_linear / max_Sxx + 1e-12)

    zoom_start = max(0, len(valid_disp_time) - int(0.1 / dt))
    steady_time = valid_disp_time[zoom_start:]
    steady_x1l = valid_x1l[zoom_start:]
    steady_x1r = valid_x1r[zoom_start:]

    save_displacement_waveform(
        pressure_val, zeta_val, steady_time, steady_x1l, steady_x1r
    )

    # =========================
    # 描画（スペクトログラムのみ）
    # =========================
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
    mesh = ax.pcolormesh(
        t_spec, f / 1000, Sxx_db, shading='gouraud', cmap='magma',
        vmin=-80, vmax=0,
    )
    ax.set_title(
        f"Airflow spectrogram (Ps = {pressure_val} Pa, zeta = {zeta_val:g})",
        fontsize=16,
    )
    ax.set_xlabel("Time [s]", fontsize=16)
    ax.set_ylabel("Frequency [kHz]", fontsize=16)
    ax.set_ylim(0, 4)
    ax.tick_params(direction='in')

    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label('SPL [dB]', fontsize=14)

    plt.tight_layout()
    
    # 画像の保存
    save_filename = OUTPUT_DIR / f"analysis_{condition_tag(pressure_val, zeta_val)}.png"
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 画像を保存しました: {save_filename}\n")


def save_displacement_waveform(
    pressure_val, zeta_val, steady_time, steady_x1l, steady_x1r
):
    fig, ax = plt.subplots(figsize=(8, 3), dpi=150)

    ax.plot(steady_time, steady_x1l, label='Left VF', color='#0072B2', linewidth=1.2)
    ax.plot(steady_time, steady_x1r, label='Right VF', color='#D55E00', linewidth=1.2, linestyle='--')
    #ax.set_title(f"Vocal Fold Displacement (Ps = {pressure_val} Pa)", fontsize=16)
    ax.set_xlabel("Time [s]", fontsize=18)
    ax.set_ylabel("Displacement [mm]", fontsize=18)
    # ax.set_ylim(-0.85, 3.3)
    ax.tick_params(direction='in')
    ax.legend(loc='upper right', fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()

    save_filename = OUTPUT_DIR / (
        f"displacement_waveform_{condition_tag(pressure_val, zeta_val)}.png"
    )
    plt.savefig(save_filename, dpi=300)
    plt.close()
    print(f"[出力完了] 声帯変位波形を保存しました: {save_filename}")


# =========================
# メインのスイープ実行ループ
# =========================
if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CSV_DIR.mkdir(parents=True, exist_ok=True)
    classification_records = []
    for p in pressure_list:
        for zeta in zeta_list:
            print(
                f"========== Start Simulation: Ps={p} Pa, "
                f"zetaL=zetaR={zeta:g} =========="
            )
            update_param_file(PARAM_FILE, p, zeta)

            print("[実行中] C++ソルバーを計算しています...")
            result = subprocess.run(
                [EXECUTABLE, PARAM_FILE], cwd=Path(EXECUTABLE).parent,
                capture_output=True, text=True,
            )

            if result.returncode != 0:
                print(
                    f"[エラー] シミュレーションが異常終了しました "
                    f"(Ps={p} Pa, zeta={zeta:g})"
                )
                print(result.stderr)
                continue
            
            analyze_and_plot(p, zeta)
            save_flow_waveform(p, zeta)
            save_full_history(p, zeta)
            classification_records.append(classify_and_save(p, zeta))
            write_classification_summary(classification_records)
        
    print("========== 全てのスイープ解析が完了しました！ ==========")

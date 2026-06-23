#!/usr/bin/env python3
"""Compare Welch-averaged power spectra from two WAV files."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.io import wavfile


WINDOWS = ("hann", "hamming", "blackman", "none")
REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare Welch-averaged power spectra from two WAV files."
    )
    parser.add_argument(
        "wav_a",
        nargs="?",
        type=Path,
        default=REPO_ROOT / "output" / "simulation.wav",
        help="First WAV file. Default: output/test_sound.wav",
    )
    parser.add_argument(
        "wav_b",
        nargs="?",
        type=Path,
        default=REPO_ROOT / "output" / "subject2.wav",
        help="Second WAV file. Default: output/subject2.wav",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("result_compare_wav_power_spectrum.png"),
        help="Output PNG path. Default: result_compare_wav_power_spectrum.png",
    )
    parser.add_argument("--label-a", default=None, help="Legend label for wav_a.")
    parser.add_argument("--label-b", default=None, help="Legend label for wav_b.")
    parser.add_argument("--start-a", type=float, default=0.0, help="Start time for wav_a [s].")
    parser.add_argument("--end-a", type=float, default=None, help="End time for wav_a [s].")
    parser.add_argument("--start-b", type=float, default=0.0, help="Start time for wav_b [s].")
    parser.add_argument("--end-b", type=float, default=None, help="End time for wav_b [s].")
    parser.add_argument(
        "--channel-a",
        type=int,
        default=None,
        help="Channel index for wav_a. Default: average all channels.",
    )
    parser.add_argument(
        "--channel-b",
        type=int,
        default=None,
        help="Channel index for wav_b. Default: average all channels.",
    )
    parser.add_argument(
        "--nperseg",
        type=int,
        default=16384,
        help="Samples per Welch segment. Default: 16384",
    )
    parser.add_argument(
        "--noverlap",
        type=int,
        default=None,
        help="Overlapped samples between segments. Default: nperseg / 2",
    )
    parser.add_argument(
        "--window",
        choices=WINDOWS,
        default="hann",
        help="Welch window. Default: hann",
    )
    parser.add_argument(
        "--average",
        choices=("mean", "median"),
        default="mean",
        help="Welch averaging method. Default: mean",
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=3000.0,
        help="Maximum plotted frequency [Hz]. Default: 3000",
    )
    parser.add_argument(
        "--db-floor",
        type=float,
        default=-140.0,
        help="Lower dB limit for plotting. Default: -140",
    )
    parser.add_argument(
        "--absolute-db",
        action="store_true",
        help="Plot dB from raw spectrum values instead of normalizing each curve to 0 dB.",
    )
    parser.add_argument(
        "--normalize-min-freq",
        type=float,
        default=20.0,
        help="Lowest frequency used for normalized dB and peak detection. Default: 20 Hz",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV output path. Frequency grids are written as separate columns.",
    )
    return parser.parse_args()


def wav_to_float(data: np.ndarray) -> np.ndarray:
    if np.issubdtype(data.dtype, np.floating):
        return data.astype(np.float64)

    if np.issubdtype(data.dtype, np.signedinteger):
        scale = float(np.iinfo(data.dtype).max)
        return data.astype(np.float64) / scale

    if np.issubdtype(data.dtype, np.unsignedinteger):
        info = np.iinfo(data.dtype)
        midpoint = (info.max + 1) / 2.0
        return (data.astype(np.float64) - midpoint) / midpoint

    raise TypeError(f"Unsupported WAV dtype: {data.dtype}")


def select_channel(data: np.ndarray, channel: int | None) -> np.ndarray:
    if data.ndim == 1:
        if channel not in (None, 0):
            raise ValueError("Mono WAV only has channel 0.")
        return data

    if channel is None:
        return np.mean(data, axis=1)

    if channel < 0 or channel >= data.shape[1]:
        raise ValueError(f"Channel must be in range 0..{data.shape[1] - 1}.")
    return data[:, channel]


def read_segment(
    wav_path: Path,
    channel: int | None,
    start: float,
    end: float | None,
) -> tuple[int, np.ndarray]:
    sample_rate, raw_data = wavfile.read(wav_path)
    samples = select_channel(wav_to_float(raw_data), channel)

    start_index = max(0, int(round(start * sample_rate)))
    end_index = len(samples) if end is None else int(round(end * sample_rate))
    end_index = min(len(samples), end_index)
    if end_index <= start_index:
        raise ValueError(f"Selected time range is empty for {wav_path}.")

    segment = samples[start_index:end_index]
    segment = segment - np.mean(segment)
    return sample_rate, segment


def welch_power(
    samples: np.ndarray,
    sample_rate: int,
    nperseg: int,
    noverlap: int | None,
    window: str,
    average: str,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    actual_nperseg = min(nperseg, len(samples))
    if actual_nperseg < 2:
        raise ValueError("At least 2 samples are required for Welch analysis.")

    actual_noverlap = actual_nperseg // 2 if noverlap is None else noverlap
    if actual_noverlap < 0 or actual_noverlap >= actual_nperseg:
        raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg.")

    scipy_window = "boxcar" if window == "none" else window
    freq, power = signal.welch(
        samples,
        fs=sample_rate,
        window=scipy_window,
        nperseg=actual_nperseg,
        noverlap=actual_noverlap,
        detrend="constant",
        return_onesided=True,
        scaling="density",
        average=average,
    )
    return freq, power, actual_nperseg, actual_noverlap


def reference_power(freq: np.ndarray, power: np.ndarray, min_freq: float) -> float:
    mask = freq >= min_freq
    if np.any(mask):
        return float(np.max(power[mask]))
    return float(np.max(power))


def peak_index(freq: np.ndarray, power: np.ndarray, min_freq: float) -> int:
    mask = freq >= min_freq
    if np.any(mask):
        masked_indices = np.where(mask)[0]
        return int(masked_indices[np.argmax(power[mask])])
    return int(np.argmax(power))


def to_db(freq: np.ndarray, power: np.ndarray, normalize: bool, min_freq: float) -> np.ndarray:
    reference = reference_power(freq, power, min_freq) if normalize else 1.0
    return 10.0 * np.log10(power / reference + 1.0e-24)


def save_csv(path: Path, freq_a: np.ndarray, db_a: np.ndarray, freq_b: np.ndarray, db_b: np.ndarray) -> None:
    n_rows = max(len(freq_a), len(freq_b))
    table = np.full((n_rows, 4), np.nan)
    table[: len(freq_a), 0] = freq_a
    table[: len(db_a), 1] = db_a
    table[: len(freq_b), 2] = freq_b
    table[: len(db_b), 3] = db_b
    np.savetxt(
        path,
        table,
        delimiter=",",
        header="freq_a_hz,power_a_db,freq_b_hz,power_b_db",
        comments="",
    )


def main() -> None:
    args = parse_args()

    fs_a, samples_a = read_segment(args.wav_a, args.channel_a, args.start_a, args.end_a)
    fs_b, samples_b = read_segment(args.wav_b, args.channel_b, args.start_b, args.end_b)

    freq_a, power_a, nperseg_a, noverlap_a = welch_power(
        samples_a, fs_a, args.nperseg, args.noverlap, args.window, args.average
    )
    freq_b, power_b, nperseg_b, noverlap_b = welch_power(
        samples_b, fs_b, args.nperseg, args.noverlap, args.window, args.average
    )

    normalize = not args.absolute_db
    db_a = to_db(freq_a, power_a, normalize, args.normalize_min_freq)
    db_b = to_db(freq_b, power_b, normalize, args.normalize_min_freq)

    label_a = args.label_a or args.wav_a.stem
    label_b = args.label_b or args.wav_b.stem

    fig, ax = plt.subplots(figsize=(8, 4), dpi=100)
    ax.plot(freq_a, db_a, color="cornflowerblue", linewidth=1.5, label=label_a)
    ax.plot(freq_b, db_b, color="firebrick", linewidth=1.5, linestyle="--", label=label_b)

    ax.set_xlabel("Frequency [Hz]", fontsize=20)
    ax.set_ylabel("Power [dB]" if args.absolute_db else "PSD [dB ]", fontsize=20)
    ax.set_xlim(0, args.max_freq)
    ax.set_ylim(-100, 5)
    ax.tick_params(direction="in", labelsize=16, top=True, right=True)
    ax.grid(which="both", linestyle="--", alpha=0.5)
    ax.legend(frameon=False, fontsize=20)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)

    if args.csv is not None:
        save_csv(args.csv, freq_a, db_a, freq_b, db_b)

    peak_a = peak_index(freq_a, power_a, args.normalize_min_freq)
    peak_b = peak_index(freq_b, power_b, args.normalize_min_freq)
    print(f"{label_a}: fs={fs_a} Hz, samples={len(samples_a)}, nperseg={nperseg_a}, noverlap={noverlap_a}")
    print(f"{label_a}: peak_frequency={freq_a[peak_a]:.2f} Hz")
    print(f"{label_b}: fs={fs_b} Hz, samples={len(samples_b)}, nperseg={nperseg_b}, noverlap={noverlap_b}")
    print(f"{label_b}: peak_frequency={freq_b[peak_b]:.2f} Hz")
    print(f"wrote: {args.output}")
    if args.csv is not None:
        print(f"wrote: {args.csv}")


if __name__ == "__main__":
    main()

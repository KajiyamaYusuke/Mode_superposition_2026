#!/usr/bin/env python3
"""Read a WAV file and output its single-sided power spectrum.

Examples:
    python tools/wav_power_spectrum.py input.wav
    python tools/wav_power_spectrum.py input.wav -o result_power.png --csv result_power.csv
    python tools/wav_power_spectrum.py input.wav --start 0.05 --end 0.4 --max-freq 6000
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile


WINDOWS = {
    "hann": np.hanning,
    "hamming": np.hamming,
    "blackman": np.blackman,
    "none": lambda n: np.ones(n),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read a WAV file and output a power spectrum plot and CSV."
    )
    parser.add_argument("wav_file", type=Path, help="Input WAV file.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("result_wav_power_spectrum.png"),
        help="Output PNG path. Default: result_wav_power_spectrum.png",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="Optional CSV output path for frequency,power,power_db.",
    )
    parser.add_argument(
        "--channel",
        type=int,
        default=None,
        help="Channel index for multi-channel WAV. Default: average all channels.",
    )
    parser.add_argument(
        "--start",
        type=float,
        default=0.0,
        help="Start time in seconds. Default: 0.0",
    )
    parser.add_argument(
        "--end",
        type=float,
        default=None,
        help="End time in seconds. Default: end of file.",
    )
    parser.add_argument(
        "--window",
        choices=sorted(WINDOWS),
        default="hann",
        help="FFT window. Default: hann",
    )
    parser.add_argument(
        "--max-freq",
        type=float,
        default=None,
        help="Maximum plotted frequency in Hz. Default: Nyquist frequency.",
    )
    parser.add_argument(
        "--db-floor",
        type=float,
        default=-140.0,
        help="Lower dB limit for plotting. Default: -140",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Plot linear power instead of normalized dB.",
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


def compute_power_spectrum(
    samples: np.ndarray, sample_rate: int, window_name: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    samples = samples - np.mean(samples)
    window = WINDOWS[window_name](len(samples))
    windowed = samples * window

    spectrum = np.fft.rfft(windowed)
    freq = np.fft.rfftfreq(len(windowed), d=1.0 / sample_rate)

    # Single-sided power spectrum with coherent window-power correction.
    window_power = np.sum(window**2)
    power = (np.abs(spectrum) ** 2) / window_power
    if len(power) > 2:
        power[1:-1] *= 2.0

    max_power = np.max(power)
    power_db = 10.0 * np.log10(power / max_power + 1.0e-24)
    return freq, power, power_db


def main() -> None:
    args = parse_args()

    sample_rate, raw_data = wavfile.read(args.wav_file)
    samples = select_channel(wav_to_float(raw_data), args.channel)

    start_index = max(0, int(round(args.start * sample_rate)))
    end_index = len(samples) if args.end is None else int(round(args.end * sample_rate))
    end_index = min(len(samples), end_index)
    if end_index <= start_index:
        raise ValueError("Selected time range is empty.")

    samples = samples[start_index:end_index]
    freq, power, power_db = compute_power_spectrum(samples, sample_rate, args.window)

    if args.csv is not None:
        csv_data = np.column_stack((freq, power, power_db))
        np.savetxt(
            args.csv,
            csv_data,
            delimiter=",",
            header="frequency_hz,power,power_db_re_max",
            comments="",
        )

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    y = power if args.linear else power_db
    ax.plot(freq, y, color="cornflowerblue", alpha=0.85)

    ax.set_xlabel("Frequency [Hz]", fontsize=16)
    ax.set_ylabel("Power" if args.linear else "Power [dB re max]", fontsize=16)
    ax.tick_params(direction="in", labelsize=12, top=True, right=True)
    ax.grid(which="both", linestyle="--", alpha=0.5)

    x_max = sample_rate / 2.0 if args.max_freq is None else args.max_freq
    ax.set_xlim(0, x_max)
    if not args.linear:
        ax.set_ylim(args.db_floor, 5)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)

    peak_index = int(np.argmax(power[1:]) + 1) if len(power) > 1 else 0
    print(f"sample_rate: {sample_rate} Hz")
    print(f"samples: {len(samples)}")
    print(f"peak_frequency: {freq[peak_index]:.2f} Hz")
    print(f"wrote: {args.output}")
    if args.csv is not None:
        print(f"wrote: {args.csv}")


if __name__ == "__main__":
    main()

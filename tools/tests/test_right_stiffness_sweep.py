#!/usr/bin/env python3

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from run_right_stiffness_sweep import make_parameter_text  # noqa: E402


def test_generated_parameter_uses_absolute_inputs_and_replaces_right_model(tmp_path):
    base_dir = tmp_path / "input"
    text = "\n".join([
        "leftFrequencyFile = M5_test/left_freq.txt",
        "rightFrequencyFile = M5_test/old_right_freq.txt",
        "leftModeFile = M5_test/left_mode.vtu",
        "rightModeFile = M5_test/old_right_mode.vtu",
        "leftSurfaceNasFile = M5_test/surface.nas",
        "rightSurfaceNasFile = M5_test/surface.nas",
    ])
    generated, frequency, mode = make_parameter_text(text, base_dir, 8, 2)
    assert f"rightFrequencyFile = {frequency}" in generated
    assert f"rightModeFile = {mode}" in generated
    assert str(base_dir.resolve() / "M5_test/left_mode.vtu") in generated
    assert frequency.name == "M5_freq_T3_d2_b8c2.txt"
    assert mode.name == "M5_mode_T3_b8c2.vtu"

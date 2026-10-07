#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_stiffness_metrics import (  # noqa: E402
    bootstrap_ratio,
    distribution_stats,
    local_variation,
    parse_runs,
)


def test_local_jitter_definition():
    periods = np.array([1.0, 1.1, 0.9])
    expected = np.mean([0.1, 0.2]) / np.mean(periods)
    assert np.isclose(local_variation(periods), expected)


def test_distribution_uses_sample_standard_deviation():
    stats = distribution_stats(np.array([1.0, 2.0, 3.0]), "value")
    assert stats["value_mean"] == 2.0
    assert stats["value_median"] == 2.0
    assert stats["value_std"] == 1.0


def test_bootstrap_ratio_is_reproducible_and_reports_interval():
    left = np.array([2.0, 2.1, 1.9, 2.2])
    right = np.array([1.0, 1.1, 0.9, 1.0])
    first = bootstrap_ratio(left, right, samples=500, rng=np.random.default_rng(7))
    second = bootstrap_ratio(left, right, samples=500, rng=np.random.default_rng(7))
    assert first == second
    assert first["ci95_low"] <= first["ratio"] <= first["ci95_high"]
    assert first["bootstrap_std"] > 0.0


def test_default_run_targets_current_output():
    assert parse_runs([]) == [("current", Path(__file__).resolve().parents[2] / "output")]

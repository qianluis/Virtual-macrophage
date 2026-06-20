"""Test chemotaxis along a CCL2 gradient."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from experiments.chemotaxis import run


def test_chemotaxis_drifts_up_gradient():
    # Average over many seeds — single random walk is too noisy
    grad = [run(seed=s, n_ticks=200, with_gradient=True)["displacement_x"]
            for s in range(20)]
    null = [run(seed=s, n_ticks=200, with_gradient=False)["displacement_x"]
            for s in range(20)]
    g_mean = float(np.mean(grad))
    n_mean = float(np.mean(null))
    # With a +x CCL2 gradient, the cell drifts up-gradient on average.
    # Without gradient, mean displacement should be small (random walk → ~0).
    assert g_mean > 5.0, f"gradient drift weak: {g_mean}"
    assert abs(n_mean) < g_mean / 2, f"null drift not zero-like: {n_mean} vs {g_mean}"

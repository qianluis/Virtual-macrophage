"""Test the LPS→M1 golden-standard behavior.

Asserts the deterministic core matches the published expectations for LPS
exposure (Mosser & Edwards 2008; O'Neill et al. 2016):
  * NF-κB activates
  * Polarization → M1 (score > 0.3)
  * Metabolism → glycolysis
  * Morphology → spread
  * TNF-α / IL-6 / IL-1β secreted; IL-10 / TGF-β NOT secreted
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def test_lps_drives_m1_polarization():
    world = World(grid_shape=(3, 3, 1))
    world.add_ligand("lps", initial=100.0)
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])

    for _ in range(60):
        cell.tick(world)

    s = cell.summary()
    # NF-κB should be activated
    assert s["tfs"]["nfkb"] > 0.5, f"NF-κB not activated: {s['tfs']['nfkb']}"
    # IRF3 (TRIF arm) also engaged
    assert s["tfs"]["irf3"] > 0.3, f"IRF3 not activated: {s['tfs']['irf3']}"
    # M2 TFs should not be activated
    assert s["tfs"]["stat6"] < 0.05
    assert s["tfs"]["stat3"] < 0.05
    # polarization → M1
    assert s["polarization_score"] > 0.3, f"not M1: {s['polarization_score']}"
    assert s["polarization_label"] == "M1"
    # metabolism switched to glycolysis
    assert s["metabolism"] == "glycolysis"
    # morphology → spread (no particles, so not cup)
    assert s["morphology"] == "spread"
    # M1 cytokines secreted
    cyto = s["cytokines_secreted_total"]
    assert cyto["tnf_alpha"] > 1.0, f"TNF-α not secreted: {cyto['tnf_alpha']}"
    assert cyto["il6"] > 0.5
    assert cyto["il1b"] > 0.3
    # M2 cytokines NOT secreted under LPS alone
    assert cyto["il10"] == 0.0
    assert cyto["tgfb"] == 0.0


def test_no_stimulus_stays_m0():
    world = World(grid_shape=(3, 3, 1))
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])

    for _ in range(60):
        cell.tick(world)

    s = cell.summary()
    assert s["polarization_label"] == "M0"
    assert s["metabolism"] == "neutral"
    assert s["morphology"] == "round"
    assert s["cytokines_secreted_total"]["tnf_alpha"] == 0.0
    assert s["cytokines_secreted_total"]["il10"] == 0.0

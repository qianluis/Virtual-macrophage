"""Test the IL-4→M2 golden-standard behavior."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def test_il4_drives_m2_polarization():
    world = World(grid_shape=(3, 3, 1))
    world.add_ligand("il4", initial=20.0)
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])

    for _ in range(60):
        cell.tick(world)

    s = cell.summary()
    assert s["tfs"]["stat6"] > 0.5
    assert s["tfs"]["nfkb"] < 0.05
    assert s["tfs"]["stat1"] < 0.05
    assert s["polarization_score"] < -0.3
    assert s["polarization_label"] == "M2"
    assert s["metabolism"] == "oxphos"
    assert s["morphology"] == "spread"
    cyto = s["cytokines_secreted_total"]
    assert cyto["il10"] > 0.5
    assert cyto["tgfb"] > 0.3
    assert cyto["tnf_alpha"] == 0.0
    assert cyto["il6"] == 0.0

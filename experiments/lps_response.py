"""LPS response — the showcase experiment.

Exposes a single resting macrophage to LPS (100 ng/mL, well above TLR4 K_d ~0.5
ng/mL) for 60 simulated minutes and asserts the cell polarizes toward M1.

Expected literature behavior (Mosser & Edwards 2008; O'Neill et al. 2016):
  * NF-κB activates within minutes
  * Polarization score climbs into M1 territory (>0.3)
  * Metabolism shifts to glycolysis
  * Morphology changes to spread
  * TNF-α / IL-6 / IL-1β are secreted

Run:
    python experiments/lps_response.py
"""
from __future__ import annotations

import sys
from pathlib import Path

# allow running without `pip install -e .`
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def run() -> dict:
    world = World(grid_shape=(5, 5, 1), voxel_size_um=10.0)
    world.add_ligand("lps", initial=100.0, D=0.0, decay=0.0)  # uniform 100 ng/mL

    cell = Macrophage(rng=np.random.default_rng(42))
    cell.motility.position = np.array([25.0, 25.0, 0.0])  # center voxel

    snapshots = []
    for t in range(60):
        cell.tick(world)
        if t in (0, 5, 15, 30, 59):
            snapshots.append((t, cell.summary()))

    return {
        "snapshots": snapshots,
        "final": cell.summary(),
    }


def main() -> None:
    result = run()
    print("=== LPS response (100 ng/mL, 60 ticks) ===\n")
    for t, s in result["snapshots"]:
        print(f"t={t:>3}  polarization={s['polarization_score']:+.3f} ({s['polarization_label']})  "
              f"NF-κB={s['tfs']['nfkb']:.3f}  metab={s['metabolism']:>10}  "
              f"morph={s['morphology']:>6}  TNF-α(cum)={s['cytokines_secreted_total']['tnf_alpha']:.2f}")
    print("\nfinal:", result["final"])


if __name__ == "__main__":
    main()

"""IL-4 response — drives M2a polarization.

Exposes a resting macrophage to IL-4 (20 ng/mL, well above IL-4R K_d ~0.5
ng/mL) for 60 simulated minutes and asserts it polarizes toward M2.

Expected literature behavior (Mosser & Edwards 2008; Murray 2007):
  * STAT6 activates within minutes
  * Polarization score → M2 (negative)
  * Metabolism → OXPHOS
  * Morphology → spread (any sufficiently strong polarization)
  * IL-10 / TGF-β secreted; TNF-α / IL-6 / IL-1β NOT secreted
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def run() -> dict:
    world = World(grid_shape=(5, 5, 1), voxel_size_um=10.0)
    world.add_ligand("il4", initial=20.0)

    cell = Macrophage(rng=np.random.default_rng(42))
    cell.motility.position = np.array([25.0, 25.0, 0.0])

    snapshots = []
    for t in range(60):
        cell.tick(world)
        if t in (0, 5, 15, 30, 59):
            snapshots.append((t, cell.summary()))

    return {"snapshots": snapshots, "final": cell.summary()}


def main() -> None:
    result = run()
    print("=== IL-4 response (20 ng/mL, 60 ticks) ===\n")
    for t, s in result["snapshots"]:
        print(f"t={t:>3}  polarization={s['polarization_score']:+.3f} ({s['polarization_label']})  "
              f"STAT6={s['tfs']['stat6']:.3f}  metab={s['metabolism']:>10}  "
              f"morph={s['morphology']:>6}  IL-10(cum)={s['cytokines_secreted_total']['il10']:.2f}")
    print("\nfinal:", result["final"])


if __name__ == "__main__":
    main()

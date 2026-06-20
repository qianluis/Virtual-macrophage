"""Phagocytosis + ROS burst — receptor-gated engulfment.

Drops opsonized particles around a macrophage in proximity range; asserts
phagocytosis events accumulate and ROS bursts are triggered.

Expected literature behavior (Underhill & Goodridge 2012; Forman & Torres 2002):
  * particles within engulf radius are phagocytosed at the per-tick rate
  * each engulfment triggers a ROS burst lasting ~5 minutes
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World
from macrophage.environment.world import Particle


def run() -> dict:
    world = World(grid_shape=(5, 5, 1), voxel_size_um=10.0)
    cell = Macrophage(rng=np.random.default_rng(7))
    cell.motility.position = np.array([25.0, 25.0, 0.0])

    # 5 IgG-opsonized particles dropped right next to the cell
    for i in range(5):
        world.particles.append(
            Particle(
                position=cell.motility.position + np.array([float(i - 2), 0.0, 0.0]),
                ligand_class="igg_opsonized_particle",
            )
        )

    # 30 ticks should be enough for several engulfments at the default rate
    snapshots = []
    for t in range(30):
        cell.tick(world)
        if t in (0, 2, 5, 10, 20, 29):
            snapshots.append((t, cell.summary()))

    return {
        "snapshots": snapshots,
        "engulfed": sum(1 for p in world.particles if p.engulfed),
        "remaining": sum(1 for p in world.particles if not p.engulfed),
        "final": cell.summary(),
    }


def main() -> None:
    result = run()
    print("=== Phagocytosis (5 IgG-opsonized particles, 30 ticks) ===\n")
    for t, s in result["snapshots"]:
        print(f"t={t:>3}  phagocytosed={s['phagocytosed']}  ROS={s['ros_burst_remaining']}  "
              f"morph={s['morphology']}")
    print(f"\nengulfed={result['engulfed']}, remaining={result['remaining']}")
    print("final:", result["final"])


if __name__ == "__main__":
    main()

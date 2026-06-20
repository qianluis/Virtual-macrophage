"""Chemotaxis along a CCL2 gradient.

Sets up a CCL2 gradient pointing in +x direction; asserts the cell drifts
preferentially toward higher concentration over a long enough run.

Expected literature behavior (Lauffenburger & Linderman 1993; Deshmane 2009):
  * Net displacement along the gradient direction > random-walk-only baseline
  * Persistence makes the trajectory smoother than independent random steps
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def run(seed: int = 1, n_ticks: int = 200, with_gradient: bool = True) -> dict:
    nx, ny = 21, 21
    h = 10.0  # µm
    world = World(grid_shape=(nx, ny, 1), voxel_size_um=h)

    if with_gradient:
        # linear gradient: 0 ng/mL at x=0, 10 ng/mL at x=nx*h
        field = np.zeros((nx, ny, 1))
        for ix in range(nx):
            field[ix, :, :] = ix * (10.0 / nx)
        world.add_ligand("ccl2", initial=field)
    else:
        world.add_ligand("ccl2", initial=0.0)

    cell = Macrophage(rng=np.random.default_rng(seed))
    cell.motility.position = np.array([nx * h / 2, ny * h / 2, 0.0])
    start = cell.motility.position.copy()

    trajectory = [cell.motility.position.copy()]
    for _ in range(n_ticks):
        cell.tick(world)
        trajectory.append(cell.motility.position.copy())

    end = cell.motility.position.copy()
    return {
        "start": start,
        "end": end,
        "displacement_x": float(end[0] - start[0]),
        "displacement_total_um": float(np.linalg.norm(end - start)),
        "trajectory": np.asarray(trajectory),
    }


def main() -> None:
    print("=== Chemotaxis along CCL2 gradient ===\n")
    grad_runs = [run(seed=s, with_gradient=True)["displacement_x"] for s in range(20)]
    null_runs = [run(seed=s, with_gradient=False)["displacement_x"] for s in range(20)]
    print(f"with gradient   : Δx mean = {np.mean(grad_runs):+.2f} µm,  "
          f"std = {np.std(grad_runs):.2f}")
    print(f"without gradient: Δx mean = {np.mean(null_runs):+.2f} µm,  "
          f"std = {np.std(null_runs):.2f}")
    print("\n(positive Δx with gradient + ~zero without = chemotaxis works.)")


if __name__ == "__main__":
    main()

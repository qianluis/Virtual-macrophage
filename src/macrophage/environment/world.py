"""Environment — the World object that the cell perceives and acts on.

v0.1 World: a small 3D voxel grid of ligand concentrations + a list of
particles (engulfable). The grid is just a dict {ligand_name: ndarray}; cells
sample local ligands at their voxel and gradients via finite differences.

Real diffusion solvers (lattice-Boltzmann etc.) come in v0.3 when populations
of cells make spatial structure interesting.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np


@dataclass
class Particle:
    """An engulfable particle (e.g. opsonized bead, fungal cell wall fragment).

    `ligand_class` flags receptor recognition: "igg_opsonized_particle" → FcγR;
    "beta_glucan" → Dectin-1; "oxidized_ldl" → CD36; "generic" → mannose-R
    (not in v0.1).
    """
    position: np.ndarray
    ligand_class: str
    radius_um: float = 0.5
    engulfed: bool = False


@dataclass
class World:
    """A simple 3D voxel world.

    grid_shape: (nx, ny, nz)
    voxel_size_um: edge length of one voxel (default 10 µm — single-cell scale)
    """
    grid_shape: tuple[int, int, int] = (20, 20, 1)
    voxel_size_um: float = 10.0
    ligands: dict[str, np.ndarray] = field(default_factory=dict)
    particles: list[Particle] = field(default_factory=list)
    diffusion_coef: dict[str, float] = field(default_factory=dict)
    decay_per_tick: dict[str, float] = field(default_factory=dict)

    def add_ligand(
        self,
        name: str,
        initial: np.ndarray | float = 0.0,
        D: float = 0.0,
        decay: float = 0.0,
    ) -> None:
        if isinstance(initial, (int, float)):
            field_arr = np.full(self.grid_shape, float(initial))
        else:
            field_arr = np.asarray(initial, dtype=float)
            assert field_arr.shape == self.grid_shape, "ligand field shape mismatch"
        self.ligands[name] = field_arr
        self.diffusion_coef[name] = D
        self.decay_per_tick[name] = decay

    def voxel_of(self, position_um: np.ndarray) -> tuple[int, int, int]:
        ix = int(np.clip(position_um[0] / self.voxel_size_um, 0, self.grid_shape[0] - 1))
        iy = int(np.clip(position_um[1] / self.voxel_size_um, 0, self.grid_shape[1] - 1))
        iz = int(np.clip(position_um[2] / self.voxel_size_um, 0, self.grid_shape[2] - 1))
        return ix, iy, iz

    def local_concentration(self, name: str, position_um: np.ndarray) -> float:
        if name not in self.ligands:
            return 0.0
        ix, iy, iz = self.voxel_of(position_um)
        return float(self.ligands[name][ix, iy, iz])

    def gradient(self, name: str, position_um: np.ndarray) -> np.ndarray:
        """Central-difference gradient at the cell's voxel; returns ng/mL per µm."""
        if name not in self.ligands:
            return np.zeros(3)
        f = self.ligands[name]
        ix, iy, iz = self.voxel_of(position_um)
        h = self.voxel_size_um
        nx, ny, nz = self.grid_shape

        def diff(axis: int, i: int, n: int) -> float:
            if n == 1:
                return 0.0
            ip, im = min(i + 1, n - 1), max(i - 1, 0)
            sl_p = [ix, iy, iz]; sl_p[axis] = ip
            sl_m = [ix, iy, iz]; sl_m[axis] = im
            return (f[tuple(sl_p)] - f[tuple(sl_m)]) / (2 * h)

        return np.array([diff(0, ix, nx), diff(1, iy, ny), diff(2, iz, nz)])

    def deposit(self, name: str, position_um: np.ndarray, amount: float) -> None:
        if name not in self.ligands:
            self.add_ligand(name, 0.0)
        ix, iy, iz = self.voxel_of(position_um)
        self.ligands[name][ix, iy, iz] += amount

    def particles_within(
        self, position_um: np.ndarray, radius_um: float
    ) -> Iterable[Particle]:
        for p in self.particles:
            if p.engulfed:
                continue
            d = np.linalg.norm(p.position - position_um)
            if d <= radius_um:
                yield p

    def diffuse(self, dt: float = 1.0) -> None:
        """Crude isotropic diffusion + decay. v0.1: simple Laplacian step.
        Stable as long as D*dt/h² < 1/6."""
        h = self.voxel_size_um
        for name, f in list(self.ligands.items()):
            D = self.diffusion_coef.get(name, 0.0)
            decay = self.decay_per_tick.get(name, 0.0)
            if D > 0:
                lap = (
                    -6.0 * f
                    + np.roll(f, 1, axis=0) + np.roll(f, -1, axis=0)
                    + np.roll(f, 1, axis=1) + np.roll(f, -1, axis=1)
                    + (np.roll(f, 1, axis=2) + np.roll(f, -1, axis=2)
                       if f.shape[2] > 1 else 0.0)
                ) / (h * h)
                f = f + D * lap * dt
            if decay > 0:
                f = f * (1.0 - decay) ** dt
            np.maximum(f, 0.0, out=f)
            self.ligands[name] = f

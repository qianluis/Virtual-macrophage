"""Test phagocytosis + ROS burst."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World
from macrophage.environment.world import Particle


def test_opsonized_particles_get_engulfed():
    world = World(grid_shape=(3, 3, 1))
    cell = Macrophage(rng=np.random.default_rng(7))
    cell.motility.position = np.array([15.0, 15.0, 0.0])

    for i in range(5):
        world.particles.append(
            Particle(
                position=cell.motility.position + np.array([float(i - 2), 0.0, 0.0]),
                ligand_class="igg_opsonized_particle",
            )
        )

    for _ in range(50):
        cell.tick(world)

    engulfed = sum(1 for p in world.particles if p.engulfed)
    assert engulfed >= 3, f"too few engulfments: {engulfed}/5"
    assert cell.summary()["phagocytosed"] == engulfed


def test_phagocytosis_triggers_ros_burst():
    world = World(grid_shape=(3, 3, 1))
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])
    world.particles.append(
        Particle(position=cell.motility.position.copy(),
                 ligand_class="igg_opsonized_particle")
    )

    saw_ros = False
    for _ in range(20):
        cell.tick(world)
        if cell.effectors.ros_burst_remaining > 0:
            saw_ros = True
    assert saw_ros, "ROS burst never triggered after phagocytosis"


def test_unrecognized_ligand_class_not_engulfed():
    world = World(grid_shape=(3, 3, 1))
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])
    world.particles.append(
        Particle(position=cell.motility.position.copy(),
                 ligand_class="unknown_alien_polymer")  # not in receptor table
    )

    for _ in range(30):
        cell.tick(world)
    assert all(not p.engulfed for p in world.particles)
    assert cell.summary()["phagocytosed"] == 0

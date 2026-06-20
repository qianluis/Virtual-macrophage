"""Motility — biased random walk with chemotactic drift.

Drift is the gradient of log(L) of a chemoattractant ligand (here CCL2 only),
scaled by χ. Random component is unbiased Gaussian. Persistence couples each
step to the previous one's direction (reduces zig-zag).

Ref: Lauffenburger & Linderman (1993).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import params


@dataclass
class MotilityState:
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))


def step(
    state: MotilityState,
    ccl2_gradient: np.ndarray,
    ccl2_local: float,
    rng: np.random.Generator,
    dt: float = 1.0,
) -> None:
    """Update position by one tick.

    `ccl2_gradient` is ∇L of CCL2 at the cell's current position (vector). If
    CCL2 is locally absent, drift is zero (no chemotaxis). The random walk
    proceeds regardless.
    """
    speed = params.get("motility", "random_speed_um_per_min")
    chi = params.get("motility", "chemotactic_coefficient")
    persist = params.get("motility", "persistence_correlation")
    kd_ccl2 = params.get("receptors", "ccr2", "kd_ng_per_ml")

    # chemotactic drift: χ * ∇log(L) / (1 + L/Kd)
    # At low L, response is roughly linear in ∇L; at saturation, it tapers.
    if ccl2_local > 0:
        drift = chi * ccl2_gradient / max(ccl2_local, 1e-9) / (1.0 + ccl2_local / kd_ccl2)
    else:
        drift = np.zeros(3)

    # random component (Gaussian, scaled to speed)
    rand = rng.standard_normal(3) * speed

    # correlated update (persistence)
    new_velocity = persist * state.velocity + (1 - persist) * (drift + rand)
    state.velocity = new_velocity
    state.position = state.position + new_velocity * dt

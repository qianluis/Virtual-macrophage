"""Morphology — discrete state + simplified geometry summary.

States: "round" (resting), "spread" (activated/migrating), "cup" (phagocytosing
an engulfable particle). v0.1 keeps geometry as a single label + an effective
radius scalar; full membrane biophysics is v0.4+ work.

Ref: Adams & Hamilton (1984) Annu Rev Immunol 2:283-318.
"""
from __future__ import annotations

from dataclasses import dataclass

from .. import params


@dataclass
class Morphology:
    state: str = "round"   # "round" | "spread" | "cup"
    radius_um: float = 10.0  # effective radius; round ~10 µm, spread ~14 µm

    def update(self, polarization_score: float, particle_in_proximity: bool) -> None:
        spread_t = params.get("morphology", "spread_threshold")
        # cup state takes precedence: if a particle is engulfable nearby and the
        # cell is at least spread, form a phagocytic cup.
        if particle_in_proximity and abs(polarization_score) > spread_t:
            self.state = "cup"
            self.radius_um = 14.0
        elif abs(polarization_score) > spread_t:
            self.state = "spread"
            self.radius_um = 14.0
        else:
            self.state = "round"
            self.radius_um = 10.0

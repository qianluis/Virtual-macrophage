"""Metabolism — binary mode flag tied to polarization.

M1 → glycolysis (Warburg-like); M2 → OXPHOS / FAO. Hysteresis avoids chatter
at the polarization switch.

Ref: O'Neill & Pearce (2016) J Exp Med 213(1):15-23.
"""
from __future__ import annotations

from dataclasses import dataclass

from .. import params


@dataclass
class MetabolicState:
    mode: str = "neutral"  # "glycolysis" | "oxphos" | "neutral"

    def update(self, polarization_score: float) -> None:
        gly = params.get("metabolism", "glycolysis_threshold")
        oxp = params.get("metabolism", "oxphos_threshold")
        hys = params.get("metabolism", "hysteresis")

        if self.mode == "glycolysis":
            if polarization_score < gly - hys:
                self.mode = "neutral"
        elif self.mode == "oxphos":
            if polarization_score > oxp + hys:
                self.mode = "neutral"
        else:  # neutral
            if polarization_score > gly:
                self.mode = "glycolysis"
            elif polarization_score < oxp:
                self.mode = "oxphos"

"""Polarization — continuous M1/M2 score in [-1, +1].

Score follows the Murray-et-al-consensus continuum: +1 fully M1, -1 fully M2,
0 resting M0. Updated each tick from TF activations; relaxes toward 0 in
absence of stimulation.

Ref: Murray et al. (2014) Immunity 41(1):14-20.
"""
from __future__ import annotations

from dataclasses import dataclass

from .. import params
from ..transcription import TFPanel


@dataclass
class PolarizationState:
    score: float = 0.0  # -1 (M2) ... +1 (M1)

    def update(self, tfs: TFPanel) -> None:
        w_nfkb = params.get("polarization", "weight_nfkb")
        w_irf3 = params.get("polarization", "weight_irf3")
        w_stat1 = params.get("polarization", "weight_stat1")
        w_stat6 = params.get("polarization", "weight_stat6")
        w_stat3 = params.get("polarization", "weight_stat3")
        relax = params.get("polarization", "relax_per_tick")

        # `drive` is the instantaneous net direction the TFs push toward;
        # we approach that target with first-order kinetics rather than adding
        # it directly each tick (otherwise score saturates to ±1 in 1-2 ticks
        # and there's no temporal dynamics worth observing).
        target = (
            w_nfkb * tfs.nfkb.active
            + w_irf3 * tfs.irf3.active
            + w_stat1 * tfs.stat1.active
            + w_stat6 * tfs.stat6.active
            + w_stat3 * tfs.stat3.active
        )
        target = max(-1.0, min(1.0, target))
        # tau ~ 10 ticks so polarization has visible kinetics over a 60-min sim
        approach_rate = 0.10
        self.score += approach_rate * (target - self.score) - relax * self.score
        self.score = max(-1.0, min(1.0, self.score))

    @property
    def label(self) -> str:
        if self.score > 0.3:
            return "M1"
        if self.score < -0.3:
            return "M2"
        return "M0"

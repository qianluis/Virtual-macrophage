"""Transcription factors — activation from upstream receptor signal.

Each TF sums its upstream receptors' bound_fraction × activation_max, then
decays each tick toward baseline. Result is `active` in [0, 1]: a coarse-
grained proxy for nuclear-translocated active TF.

Refs: O'Neill, Golenbock & Bowie (2016) Nat Rev Immunol; Murray (2007)
Curr Opin Pharmacol (per-TF specificity).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from .. import params
from ..sensing import ReceptorPanel


@dataclass
class TF:
    name: str
    upstream: list[str]
    activation_max: float
    decay_per_tick: float
    active: float = 0.0

    def update(self, panel: ReceptorPanel) -> None:
        # signal = max(upstream receptor bound_fraction) — TFs respond to
        # whichever upstream is most engaged, not the sum (avoid runaway).
        signal = 0.0
        for rname in self.upstream:
            r = getattr(panel, rname, None)
            if r is not None:
                signal = max(signal, r.bound_fraction)
        target = signal * self.activation_max
        # first-order relaxation toward target with `decay` setting the
        # downward-only timescale; upward step is instantaneous (signal-driven)
        if target >= self.active:
            self.active = target
        else:
            self.active = max(target, self.active - self.decay_per_tick)


@dataclass
class TFPanel:
    nfkb: TF = field(default=None)  # type: ignore[assignment]
    irf3: TF = field(default=None)  # type: ignore[assignment]
    stat1: TF = field(default=None)  # type: ignore[assignment]
    stat6: TF = field(default=None)  # type: ignore[assignment]
    stat3: TF = field(default=None)  # type: ignore[assignment]

    @classmethod
    def from_params(cls) -> "TFPanel":
        def _make(name: str) -> TF:
            return TF(
                name=name,
                upstream=list(params.get("transcription", name, "upstream")),
                activation_max=params.get("transcription", name, "activation_max"),
                decay_per_tick=params.get("transcription", name, "decay_per_tick"),
            )
        p = cls()
        for n in ("nfkb", "irf3", "stat1", "stat6", "stat3"):
            setattr(p, n, _make(n))
        return p

    def update(self, panel: ReceptorPanel) -> None:
        for tf in (self.nfkb, self.irf3, self.stat1, self.stat6, self.stat3):
            tf.update(panel)

    def __iter__(self) -> Iterable[TF]:
        return iter((self.nfkb, self.irf3, self.stat1, self.stat6, self.stat3))

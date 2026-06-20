"""Macrophage — the agent.

One perceive → decide → act loop per tick. Composes the receptor panel, TF
panel, polarization/metabolism/morphology state, and effector outputs.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import params
from .advisor import MacrophageAdvisor, NullAdvisor
from .effectors import EffectorState, secretion_rates
from .environment.world import World
from .metabolism import MetabolicState
from .morphology import Morphology
from .motility import MotilityState, step as motility_step
from .polarization import PolarizationState
from .sensing import ReceptorPanel
from .transcription import TFPanel


@dataclass
class Macrophage:
    receptors: ReceptorPanel = field(default_factory=ReceptorPanel.from_params)
    tfs: TFPanel = field(default_factory=TFPanel.from_params)
    polarization: PolarizationState = field(default_factory=PolarizationState)
    metabolism: MetabolicState = field(default_factory=MetabolicState)
    morphology: Morphology = field(default_factory=Morphology)
    motility: MotilityState = field(default_factory=MotilityState)
    effectors: EffectorState = field(default_factory=EffectorState)
    advisor: MacrophageAdvisor = field(default_factory=NullAdvisor)
    rng: np.random.Generator = field(default_factory=lambda: np.random.default_rng(0))

    # ---- perceive --------------------------------------------------------
    def perceive(self, world: World) -> None:
        """Sample local ligand concentrations + update receptor binding."""
        local = {}
        # only sample ligands the receptor panel cares about
        for r in (self.receptors.tlr4, self.receptors.tlr2, self.receptors.dectin1,
                  self.receptors.cd36, self.receptors.ccr2, self.receptors.ifngr,
                  self.receptors.il4r, self.receptors.il10r):
            local[r.ligand] = world.local_concentration(r.ligand, self.motility.position)
        self.receptors.update_from_environment(local)

    # ---- decide ----------------------------------------------------------
    def decide(self, world: World | None = None) -> None:
        """Update internal state from receptor signals."""
        self.tfs.update(self.receptors)
        self.polarization.update(self.tfs)
        self.metabolism.update(self.polarization.score)

        # particle-in-proximity flag for morphology
        particle_in_proximity = False
        if world is not None:
            engulf_radius = params.get("effectors", "phagocytosis", "proximity_um")
            for _ in world.particles_within(self.motility.position, engulf_radius):
                particle_in_proximity = True
                break
        self.morphology.update(self.polarization.score, particle_in_proximity)

    # ---- act -------------------------------------------------------------
    def act(self, world: World) -> None:
        """Mutate the world: secrete cytokines, attempt phagocytosis, move."""
        # 1. cytokine secretion
        rates = secretion_rates(self.polarization.score)
        self.effectors.add_secretion(rates)
        for name, amount in (
            ("tnf_alpha", rates.tnf_alpha),
            ("il6", rates.il6),
            ("il1b", rates.il1b),
            ("il10", rates.il10),
            ("tgfb", rates.tgfb),
        ):
            if amount > 0:
                world.deposit(name, self.motility.position, amount)

        # 2. phagocytosis: engulf a nearby particle if a receptor recognizes it
        engulf_radius = params.get("effectors", "phagocytosis", "proximity_um")
        engulf_p = params.get("effectors", "phagocytosis", "rate_per_particle_per_tick")
        for p in world.particles_within(self.motility.position, engulf_radius):
            # receptor gating: only engulf if the corresponding receptor is bound
            if not self._has_engulf_receptor_for(p.ligand_class):
                continue
            if self.rng.random() < engulf_p:
                p.engulfed = True
                self.effectors.record_phagocytosis()

        # 3. ROS decay
        self.effectors.tick_ros()

        # 4. motility — chemotaxis along CCL2 gradient
        ccl2_grad = world.gradient("ccl2", self.motility.position)
        ccl2_local = world.local_concentration("ccl2", self.motility.position)
        motility_step(self.motility, ccl2_grad, ccl2_local, self.rng)

    # ---- helpers ---------------------------------------------------------
    def _has_engulf_receptor_for(self, ligand_class: str) -> bool:
        """Receptor-gating: phagocytosis requires the matching receptor to be
        engaged. v0.1 is conservative — needs explicit ligand-class match.

        Ref: Underhill & Goodridge (2012) Nat Rev Immunol 12(7):492-502.
        """
        if ligand_class == "igg_opsonized_particle":
            # FcγR not modeled as a soluble-ligand receptor; treat opsonized
            # particles as always engageable when in proximity (the receptor
            # is multivalent and avidity-driven, not concentration-driven).
            return True
        if ligand_class == "beta_glucan":
            return self.receptors.dectin1.bound_fraction > 0.05
        if ligand_class == "oxidized_ldl":
            return self.receptors.cd36.bound_fraction > 0.05
        return False  # unknown ligand class → no engulfment

    # ---- one full tick ---------------------------------------------------
    def tick(self, world: World) -> None:
        self.perceive(world)
        self.decide(world)
        self.act(world)

    # ---- summary --------------------------------------------------------
    def summary(self) -> dict[str, Any]:
        return {
            "polarization_score": round(self.polarization.score, 3),
            "polarization_label": self.polarization.label,
            "metabolism": self.metabolism.mode,
            "morphology": self.morphology.state,
            "tfs": {
                "nfkb": round(self.tfs.nfkb.active, 3),
                "stat1": round(self.tfs.stat1.active, 3),
                "stat6": round(self.tfs.stat6.active, 3),
                "stat3": round(self.tfs.stat3.active, 3),
                "irf3": round(self.tfs.irf3.active, 3),
            },
            "phagocytosed": self.effectors.phagocytosed_count,
            "ros_burst_remaining": self.effectors.ros_burst_remaining,
            "cytokines_secreted_total": {
                "tnf_alpha": round(self.effectors.cytokines_secreted_total.tnf_alpha, 3),
                "il6": round(self.effectors.cytokines_secreted_total.il6, 3),
                "il1b": round(self.effectors.cytokines_secreted_total.il1b, 3),
                "il10": round(self.effectors.cytokines_secreted_total.il10, 3),
                "tgfb": round(self.effectors.cytokines_secreted_total.tgfb, 3),
            },
            "position_um": [round(x, 2) for x in self.motility.position.tolist()],
        }

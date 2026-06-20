"""Effector outputs — cytokines, ROS, phagocytosis.

Cytokine secretion rates scale with polarization (M1 → TNF-α/IL-6/IL-1β,
M2 → IL-10/TGF-β). ROS bursts on a phagocytosis event.

Refs: Mosser & Edwards (2008); Arango Duque & Descoteaux (2014); Forman &
Torres (2002).
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .. import params


@dataclass
class CytokineOutput:
    """Per-tick cytokine secretion rates the cell will deposit into its
    environment. Units: ng/mL per tick into the local voxel."""
    tnf_alpha: float = 0.0
    il6: float = 0.0
    il1b: float = 0.0
    il10: float = 0.0
    tgfb: float = 0.0


def secretion_rates(polarization_score: float) -> CytokineOutput:
    """Map polarization score to cytokine rates. Linear in score magnitude;
    M1 cytokines scale with positive score, M2 with negative score.

    Ref: Arango Duque & Descoteaux (2014) Front Immunol 5:491.
    """
    s = max(-1.0, min(1.0, polarization_score))
    m1 = max(0.0, s)
    m2 = max(0.0, -s)
    out = CytokineOutput()
    out.tnf_alpha = m1 * params.get("effectors", "cytokines", "tnf_alpha_m1_max")
    out.il6 = m1 * params.get("effectors", "cytokines", "il6_m1_max")
    out.il1b = m1 * params.get("effectors", "cytokines", "il1b_m1_max")
    out.il10 = m2 * params.get("effectors", "cytokines", "il10_m2_max")
    out.tgfb = m2 * params.get("effectors", "cytokines", "tgfb_m2_max")
    return out


@dataclass
class EffectorState:
    """Mutable state used by effector actions."""
    ros_burst_remaining: int = 0
    phagocytosed_count: int = 0
    cytokines_secreted_total: CytokineOutput = field(default_factory=CytokineOutput)

    def tick_ros(self) -> None:
        if self.ros_burst_remaining > 0:
            self.ros_burst_remaining -= 1

    def trigger_ros_burst(self) -> None:
        # duration set on each fresh trigger; bursts overlap by reset
        self.ros_burst_remaining = max(
            self.ros_burst_remaining,
            int(params.get("effectors", "ros", "duration_ticks")),
        )

    def record_phagocytosis(self) -> None:
        self.phagocytosed_count += 1
        self.trigger_ros_burst()

    def add_secretion(self, output: CytokineOutput) -> None:
        self.cytokines_secreted_total.tnf_alpha += output.tnf_alpha
        self.cytokines_secreted_total.il6 += output.il6
        self.cytokines_secreted_total.il1b += output.il1b
        self.cytokines_secreted_total.il10 += output.il10
        self.cytokines_secreted_total.tgfb += output.tgfb

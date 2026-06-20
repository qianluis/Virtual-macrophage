"""LLM advisor hook (optional, abstract base only in v0.1).

The deterministic core handles every named stimulus the literature covers. For
**novel** stimuli the rules don't recognize — a synthetic ligand, an unusual
cytokine cocktail, an in silico drug — the user can plug in a `MacrophageAdvisor`
that proposes a response. The core then sanity-checks the proposal against the
same parameter ranges before applying.

This mirrors the nature-skills "deterministic gates first, LLM advisory second"
philosophy. v0.1 ships the ABC; v0.2+ ships an example Anthropic-API impl.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass
class AdvisorRequest:
    """Context the advisor sees when called."""
    novel_ligand: str
    local_concentration: float        # ng/mL or particle count
    current_polarization: float        # [-1, +1]
    current_morphology: str
    nearby_particles: int


@dataclass
class AdvisorResponse:
    """The advisor's proposed effect on the cell. The deterministic core
    *applies* this only after gating against `parameters.toml` ranges."""
    nudge_polarization: float = 0.0   # additive, [-0.3, 0.3] hard-clamped by core
    treat_as_known_ligand: str | None = None  # remap novel → known if defensible
    rationale: str = ""                # logged for audit


class MacrophageAdvisor(ABC):
    """Plug-in point. Subclass and implement `propose()`."""

    @abstractmethod
    def propose(self, req: AdvisorRequest) -> AdvisorResponse: ...


class NullAdvisor(MacrophageAdvisor):
    """Default: no advice. Novel stimulus → no signal (resting baseline)."""
    def propose(self, req: AdvisorRequest) -> AdvisorResponse:  # noqa: ARG002
        return AdvisorResponse(rationale="null advisor: novel stimulus ignored")

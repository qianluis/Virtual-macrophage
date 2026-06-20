"""Sensing — receptor-binding module.

Each receptor binds one ligand class and reports a `bound_fraction` in [0,1] via
a Hill function of local ligand concentration:

    bound = L^n / (Kd^n + L^n)

Receptor parameters live in data/parameters.toml under [receptors.<name>].

Ref: Janeway's Immunobiology Ch.3 (receptor recognition); per-receptor refs in
parameters.toml.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .. import params


def hill(ligand_conc: float, kd: float, n: float) -> float:
    """Hill function in [0, 1]. Ligand and Kd in same units."""
    if ligand_conc <= 0:
        return 0.0
    L = ligand_conc
    return (L ** n) / (kd ** n + L ** n)


@dataclass
class Receptor:
    name: str
    ligand: str
    kd: float
    hill_n: float
    bound_fraction: float = 0.0

    def update(self, local_ligand_conc: float) -> None:
        self.bound_fraction = hill(local_ligand_conc, self.kd, self.hill_n)


@dataclass
class ReceptorPanel:
    """The cell's six receptors used in v0.1.

    Adding a receptor: add a `[receptors.<name>]` block to parameters.toml,
    add a field here, wire it in `from_params()`. No invented Kd values.
    """
    tlr4: Receptor = field(default=None)  # type: ignore[assignment]
    tlr2: Receptor = field(default=None)  # type: ignore[assignment]
    dectin1: Receptor = field(default=None)  # type: ignore[assignment]
    cd36: Receptor = field(default=None)  # type: ignore[assignment]
    ccr2: Receptor = field(default=None)  # type: ignore[assignment]
    ifngr: Receptor = field(default=None)  # type: ignore[assignment]
    il4r: Receptor = field(default=None)  # type: ignore[assignment]
    il10r: Receptor = field(default=None)  # type: ignore[assignment]

    @classmethod
    def from_params(cls) -> "ReceptorPanel":
        def _make(name: str) -> Receptor:
            kd_key = "kd_ng_per_ml"  # all soluble ligands use ng/mL in v0.1
            return Receptor(
                name=name,
                ligand=params.get("receptors", name, "ligand"),
                kd=params.get("receptors", name, kd_key),
                hill_n=params.get("receptors", name, "hill_n"),
            )
        panel = cls()
        for name in ("tlr4", "tlr2", "dectin1", "cd36", "ccr2",
                     "ifngr", "il4r", "il10r"):
            setattr(panel, name, _make(name))
        return panel

    def update_from_environment(self, ligand_conc: dict[str, float]) -> None:
        """ligand_conc keys are ligand names from parameters.toml; absent →
        treated as 0 (no signal). An invented ligand the panel doesn't know
        is silently ignored — the caller can route unknowns to the LLM
        advisor."""
        for r in (self.tlr4, self.tlr2, self.dectin1, self.cd36,
                  self.ccr2, self.ifngr, self.il4r, self.il10r):
            r.update(ligand_conc.get(r.ligand, 0.0))

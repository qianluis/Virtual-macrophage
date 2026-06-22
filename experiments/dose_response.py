"""Dose-response — LPS concentration vs M1 polarization at steady state.

Sweeps LPS from 0.001 to 1000 ng/mL and records the steady-state polarization
score + NF-κB activity. The literature expectation (Akashi-Takamura 2008;
Mosser & Edwards 2008) is a sigmoidal dose-response: negligible below ~0.1
ng/mL, half-maximal around the TLR4 K_d (~0.5 ng/mL), saturated by ~10 ng/mL.

This is the *quantitative* companion to lps_response.py (which is the
qualitative golden standard at one saturating dose).

Run:
    python experiments/dose_response.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from macrophage import Macrophage, World


def steady_state(lps_ng_per_ml: float, n_ticks: int = 120, seed: int = 0) -> dict:
    world = World(grid_shape=(3, 3, 1))
    world.add_ligand("lps", initial=lps_ng_per_ml)
    cell = Macrophage(rng=np.random.default_rng(seed))
    cell.motility.position = np.array([15.0, 15.0, 0.0])
    for _ in range(n_ticks):
        cell.tick(world)
    s = cell.summary()
    return {
        "lps": lps_ng_per_ml,
        "polarization": s["polarization_score"],
        "nfkb": s["tfs"]["nfkb"],
        "tnf_alpha_total": s["cytokines_secreted_total"]["tnf_alpha"],
    }


def run() -> list[dict]:
    doses = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]
    return [steady_state(d) for d in doses]


def main() -> None:
    print("=== LPS dose-response (steady state at 120 ticks) ===\n")
    print(f"{'LPS (ng/mL)':>12}  {'polarization':>13}  {'NF-κB':>7}  {'TNF-α (cum)':>12}")
    for r in run():
        print(f"{r['lps']:>12.3g}  {r['polarization']:>+13.3f}  {r['nfkb']:>7.3f}  {r['tnf_alpha_total']:>12.2f}")


if __name__ == "__main__":
    main()

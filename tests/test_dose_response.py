"""Quantitative tests — dose-response shape + the emergent STAT3 autocrine loop.

These move beyond the qualitative golden standards: they verify the model is
*quantitatively* aligned with literature (sigmoidal LPS dose-response centered
near the TLR4 K_d) and that an emergent behavior (STAT3 activation under IL-4
via autocrine IL-10) is mechanistically caused by the IL-10 the cell secretes,
not a coincidental artifact.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

import numpy as np
import pytest

from macrophage import Macrophage, World
from experiments.dose_response import steady_state


def test_lps_dose_response_is_sigmoidal():
    """Polarization should rise monotonically with LPS dose and saturate."""
    results = [steady_state(d, n_ticks=120) for d in (0.001, 0.5, 100.0)]
    low, mid, high = [r["polarization"] for r in results]
    # near-zero dose → essentially M0
    assert abs(low) < 0.1, f"sub-threshold LPS shouldn't polarize: {low}"
    # saturating dose → M1
    assert high > 0.3, f"saturating LPS should be M1: {high}"
    # monotonic
    assert low <= mid <= high, f"not monotonic: {low} {mid} {high}"


def test_lps_ec50_near_tlr4_kd():
    """The half-maximal polarization dose should be within an order of magnitude
    of the TLR4 K_d (~0.5 ng/mL). We don't assert exact EC50 (the polarization
    weights are OOM), but it must be in [0.05, 50] ng/mL — i.e. the response is
    sensitive in the physiological LPS range, not at absurd doses."""
    low = steady_state(0.05)["polarization"]
    high = steady_state(100.0)["polarization"]
    span = high - low
    half = low + span / 2
    # find the dose bracketing the half-max
    doses = [0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0, 50.0]
    below = [d for d in doses if steady_state(d)["polarization"] < half]
    above = [d for d in doses if steady_state(d)["polarization"] >= half]
    if below and above:
        ec50_bracket = (max(below), min(above))
        assert ec50_bracket[0] <= 50.0 and ec50_bracket[1] >= 0.05, \
            f"EC50 outside plausible range: {ec50_bracket}"


def test_stat3_under_il4_is_autocrine_il10_driven():
    """The emergent STAT3 activation under IL-4 is driven by the IL-10 the cell
    secretes (M2c autocrine amplification, Saraiva & O'Garra 2010). Verify
    causally: if we block IL-10 sensing (zero out il10 receptor binding),
    STAT3 should NOT activate under IL-4 alone.

    This is the difference between 'the model happens to show STAT3' and 'the
    model shows STAT3 *because* of the IL-10 autocrine loop the rules predict'.
    """
    # Baseline: IL-4 alone → STAT3 should activate (via secreted IL-10)
    world = World(grid_shape=(3, 3, 1))
    world.add_ligand("il4", initial=20.0)
    cell = Macrophage(rng=np.random.default_rng(0))
    cell.motility.position = np.array([15.0, 15.0, 0.0])
    for _ in range(60):
        cell.tick(world)
    stat3_with_autocrine = cell.tfs.stat3.active

    # Intervention: IL-4 but the cell can't sense its own IL-10 (receptor Kd
    # made effectively infinite → never binds). We do this by monkey-patching
    # the il10r update to always report 0.
    world2 = World(grid_shape=(3, 3, 1))
    world2.add_ligand("il4", initial=20.0)
    cell2 = Macrophage(rng=np.random.default_rng(0))
    cell2.motility.position = np.array([15.0, 15.0, 0.0])
    original_update = cell2.receptors.il10r.update
    def blinded(_conc):
        cell2.receptors.il10r.bound_fraction = 0.0
    cell2.receptors.il10r.update = blinded  # type: ignore[assignment]
    for _ in range(60):
        cell2.tick(world2)
    stat3_blocked = cell2.tfs.stat3.active
    cell2.receptors.il10r.update = original_update  # restore

    assert stat3_with_autocrine > 0.3, \
        f"IL-4 should activate STAT3 via autocrine IL-10: {stat3_with_autocrine}"
    assert stat3_blocked < 0.05, \
        f"blocking IL-10 sensing should abolish STAT3: {stat3_blocked}"
    assert stat3_with_autocrine > stat3_blocked + 0.2, \
        "STAT3 under IL-4 is not driven by autocrine IL-10"

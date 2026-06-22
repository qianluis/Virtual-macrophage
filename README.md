# virtual-macrophage

> A virtual human macrophage as an agent — single-cell biology rules + a perceive/decide/act loop, with an optional LLM advisory layer for novel stimuli.

**Status: v0.1 (skeleton)** — single-cell, deterministic core, ~10 ligands, M1/M2 polarization, phagocytosis, chemotaxis, cytokine secretion. **Not** a whole-genome / molecular-dynamics / multicellular tissue model.

This project takes the philosophy from [`nature-skills`](https://github.com/Yuan1z0825/nature-skills) — ground every claim in literature, surface uncertainty as data, anti-fabricate via deterministic gates — and applies it to systems biology: every module's rules trace to a cited paper, every magic-number lives in `data/parameters.toml` with a citation, and the agent's behavior under named stimuli is verified against published expectations by tests.

## Why "agent", not "ODE model"

Traditional macrophage models are systems of ODEs (continuous concentrations, no notion of "this cell, right now"). That's the right tool for population-level cytokine dynamics; it's the wrong tool for asking *"what would this individual cell do if I dropped LPS + IL-4 on it simultaneously and a pathogen showed up at the north edge of its environment?"* The agent paradigm gives:

- **State per cell** — polarization vector, transcription factor activity, metabolic mode, position, cell-cycle phase.
- **Perceive → decide → act** loop — receptors sense local ligand concentrations, internal logic updates state, the cell secretes / phagocytoses / migrates as actions on the environment.
- **Composability with LLMs** — the deterministic core handles named stimuli; an `advisor.py` ABC lets a language model propose responses to **novel** stimuli not covered by the literature rules (with the deterministic core still gating what actually happens).

## What "as human as possible" means here (v0.1 scope)

A v0.1 macrophage that the literature would call recognizably human-like along these axes:

| Axis | Human macrophage | v0.1 implementation | Reference |
|---|---|---|---|
| Morphology | round (resting) → spread/elongated (activated) → cup-shaped (phagocytosing) | discrete morphology states + simplified geometry | Adams & Hamilton 1984 |
| Sensing | TLR2/4, Dectin-1, CD36, FcγR, mannose receptor, CCR2 | 6-receptor module with literature `K_d` and ligand specificities | Janeway's Immunobiology Ch.3 |
| Polarization | M0 / M1 (LPS, IFN-γ) / M2a (IL-4) / M2c (IL-10) | continuous M1/M2 score with stimulus-driven dynamics | Mosser & Edwards 2008; Murray et al. 2014 |
| Transcription | NF-κB, STAT1, STAT3, STAT6, IRF3 | simplified Hill-function activations driving cytokine output | O'Neill et al. 2016 |
| Metabolism | M1 → glycolysis / M2 → OXPHOS | binary mode flag coupled to polarization | O'Neill & Pearce 2016 |
| Effectors | phagocytosis, ROS burst, NO, cytokines (TNF-α, IL-6, IL-1β, IL-10, TGF-β) | rate-based actions on the environment | Mosser & Edwards 2008 |
| Motility | chemotaxis along CCL2/CXCL gradients + random | biased random walk with chemotactic drift | Lauffenburger & Linderman 1993 |

Out of scope for v0.1 (deferred): full GRN, membrane biophysics, MHC-II antigen presentation kinetics, tissue-scale crowds.

## What works in v0.1

A single-cell perceive/decide/act loop with **four golden-standard experiments** verified by tests:

1. `lps_response.py` — LPS (100 ng/mL) drives M1 polarization, NF-κB activation, glycolytic shift, TNF-α/IL-6 secretion.
2. `il4_response.py` — IL-4 (20 ng/mL) drives M2a polarization, STAT6 activation, OXPHOS, IL-10/TGF-β secretion.
3. `phagocytosis.py` — opsonized particles in proximity → phagocytosis events + ROS burst.
4. `chemotaxis.py` — CCL2 gradient → biased migration toward source, with directional persistence.
5. `dose_response.py` — LPS dose-response curve, verifying sigmoidal TLR4/NF-κB behavior and EC50 plausibility.

Each has a `tests/test_*.py` that runs the experiment headlessly and asserts the expected literature behavior. See [`docs/v0.1-results.md`](docs/v0.1-results.md) for the current numbers and generated figures.

## Install

```bash
git clone https://github.com/<you>/virtual-macrophage.git
cd virtual-macrophage
pip install -e '.[dev]'
pytest                                  # all tests green
python experiments/lps_response.py      # showcase: LPS -> M1
python experiments/il4_response.py      # IL-4 -> M2a + autocrine IL-10/STAT3
python experiments/dose_response.py     # LPS dose-response
PYTHONPATH=src python -m macrophage.viz # optional figures
```

Python 3.10+. Core deps: `numpy`. Optional: `matplotlib` for the experiment plots, an LLM client for the advisor.

## Anti-fabrication discipline (borrowed from nature-skills)

- Every parameter in `data/parameters.toml` carries a `source = "..."` field with a real citation. No invented K_d / EC50 / half-life values.
- Every rule in code carries a `# Ref:` comment to the paper it implements.
- The deterministic core never generates novel biology — its behavior is bounded by the cited rules. Novel stimuli either fall back to a documented default ("if unrecognized, behave as resting M0 with a chemokine sniff") or, if the `MacrophageAdvisor` is wired, route to an LLM advisory call whose output is sanity-checked against the same parameter ranges before being applied.
- A `MISSING_PARAMETER` is a placeholder, never a guess.

See [`AGENTS.md`](AGENTS.md) for the architectural overview and [`docs/biology.md`](docs/biology.md) for the literature grounding.

## License

MIT. Cite the references in `docs/biology.md` if you publish work using this.

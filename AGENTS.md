# AGENTS.md — virtual-macrophage architecture

The cell is an agent. One **perceive → decide → act** loop per simulation tick, mutating both internal state and the environment.

```
                ┌─────────────────────────────────────────────┐
                │          environment (World)                │
                │  ligand fields, particles, other cells      │
                └────────────────┬────────────────────────────┘
                                 │ local sample
                          ┌──────▼──────┐
                          │   sensing/  │  TLR2/4, Dectin-1, CD36, FcγR,
                          │  receptors  │  mannose-R, CCR2 → bound_fraction
                          └──────┬──────┘
                                 │
                ┌────────────────▼─────────────────┐
                │       transcription (TFs)         │  NF-κB, STAT1/3/6, IRF3
                │   activation = Hill(receptor)     │  Refs: O'Neill 2016
                └────────────────┬──────────────────┘
                                 │
                ┌────────────────▼─────────────────┐
                │   polarization (M0/M1/M2 score)   │  continuous in [-1, 1]
                │   updated from TF activity        │  Refs: Murray 2014
                └────────────────┬──────────────────┘
                                 │
                ┌────────────────▼─────────────────┐
                │     metabolism (mode switch)      │  glycolysis ↔ OXPHOS
                │     follows polarization          │  Refs: O'Neill 2016
                └────────────────┬──────────────────┘
                                 │
                ┌────────────────▼─────────────────┐
                │  morphology (round/spread/cup)    │  visual + functional state
                └────────────────┬──────────────────┘
                                 │
                ┌────────────────▼─────────────────┐
                │            effectors              │  phagocytosis, ROS burst,
                │    actions on the environment     │  cytokine secretion (TNF-α,
                │                                   │  IL-6, IL-1β, IL-10, TGF-β)
                └────────────────┬──────────────────┘
                                 │
                ┌────────────────▼─────────────────┐
                │             motility              │  biased random walk along
                │  chemotaxis along ligand grads    │  CCL2 / CXCL gradients
                └───────────────────────────────────┘
```

## Module map

| Module | What it does | Cited refs |
|---|---|---|
| `morphology/` | discrete state (round/spread/cup) + simplified 3D geometry summary | Adams & Hamilton 1984 |
| `sensing/` | receptor binding (Hill), local ligand sampling | Janeway Ch.3 |
| `transcription/` | TF activation from receptor signal | O'Neill 2016 |
| `polarization/` | M1/M2 score from TF activity | Mosser & Edwards 2008; Murray 2014 |
| `metabolism/` | glycolysis/OXPHOS switch tied to polarization | O'Neill & Pearce 2016 |
| `effectors/` | phagocytosis, ROS, cytokine secretion | Mosser & Edwards 2008 |
| `motility/` | biased random walk + chemotactic drift | Lauffenburger & Linderman 1993 |
| `environment/` | World tick loop, ligand fields, particles | (this project) |
| `advisor/` | optional LLM hook for novel stimuli | (this project) |

## The deterministic core / LLM advisor split

Mirrors `nature-skills`' "structural rules first, LLM second":

- **Deterministic core** (`Macrophage.tick()`): handles every named stimulus the literature covers. Behavior is bounded by parameters in `data/parameters.toml`. Tests verify it under golden-standard experiments.
- **LLM advisor** (`advisor.MacrophageAdvisor` ABC): called only when the core encounters a stimulus it has no rule for (e.g. a synthetic ligand the user invented). The advisor proposes a response; the core gates it against the same parameter ranges before applying. The advisor never bypasses the deterministic checks.

Default `advisor=None` → unrecognized stimuli fall back to documented defaults ("treat as no signal; resting baseline + chemokine sniff").

## Tick loop

```python
for t in range(n_ticks):
    cell.perceive(world)        # sample local ligand concentrations
    cell.decide()               # update receptors → TFs → polarization → metabolism → morphology
    cell.act(world)             # secrete, phagocytose, migrate; mutate world
    world.diffuse(dt)           # ligand diffusion + decay
```

A tick is conceptually `~1 minute` of biological time — fast enough to resolve cytokine secretion (minutes) and chemotaxis (minutes), slow enough to coarse-grain receptor binding (seconds).

## Parameters

Every numeric constant lives in `data/parameters.toml`:

```toml
[receptors.tlr4]
ligand = "lps"
kd_ng_per_ml = 0.5                # Akashi-Takamura 2008
hill_n = 1.5                       # fit to dose-response in same paper
source = "Akashi-Takamura & Miyake (2008) Curr Opin Immunol 20(4):420-425"
```

No invented numbers. Where a parameter is unknown / paper-specific, mark `value = "MISSING"` and the loader raises `MISSING_PARAMETER` rather than substituting a guess.

## Anti-fabrication discipline

- Every rule in code: `# Ref: <citation>` comment.
- Every parameter: `source = "..."` in TOML.
- Novel ligand without a TOML entry → log a warning, fall back to "no signal" (or call advisor).
- Tests verify the expected qualitative behavior (M1 polarization, cytokine output direction) under named stimuli.
- A `pytest -q` green run is the equivalent of `nature-skills`' `check_all` — single proof the model still matches the cited literature.

## Roadmap

- **v0.1** (this commit): single cell, 4 golden-standard experiments, deterministic core.
- **v0.2**: parameter calibration to actual published dose-response curves (RAW264.7 / THP-1 / human MDM data).
- **v0.3**: multi-cell ABM (population dynamics, paracrine signaling, M1/M2 mixed populations, simple wound).
- **v0.4**: T-cell coupling (MHC-II antigen presentation kinetics).
- **v0.5+**: tissue-scale (vasculature, fibrosis, granuloma).

Each version states its boundary explicitly, the way the research-planner skill demands.

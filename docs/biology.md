# Biological grounding

This document lists the literature each module's rules trace to. New rules go
through this gate: if a behavior isn't in the literature (or in the user's
explicit input), it doesn't go in.

## Polarization

- **Mosser & Edwards (2008)** *Nat Rev Immunol* 8(12):958-969 — the M1/M2/regulatory framework. Source for: M1 driven by IFN-γ + TLR ligands; M2a driven by IL-4 / IL-13; M2c driven by IL-10 / glucocorticoids; M1 → pro-inflammatory cytokines + ROS/NO; M2 → anti-inflammatory + tissue repair.
- **Murray et al. (2014)** *Immunity* 41(1):14-20 — the consensus nomenclature; macrophage activation is a **continuum** not discrete buckets. Source for: continuous M1/M2 score in `polarization/state.py` rather than enums.
- **Wynn, Chawla & Pollard (2013)** *Nature* 496(7446):445-455 — tissue context shapes activation. Used as the v0.2+ rationale for adding tissue-context modifiers.

## TLR signaling and TFs

- **O'Neill, Golenbock & Bowie (2016)** *Nat Rev Immunol* 13(6):453-460 — TLR signaling cascades: MyD88 → NF-κB → pro-inflammatory cytokines; TRIF → IRF3 → type-I IFN. Source for: NF-κB, IRF3 activation rules in `transcription/`.
- **Akashi-Takamura & Miyake (2008)** *Curr Opin Immunol* 20(4):420-425 — TLR4/MD-2/CD14 LPS recognition; LPS K_d ~0.5 ng/mL range. Source for: TLR4 K_d in `parameters.toml`.
- **Brown (2006)** *Nat Rev Immunol* 6(1):33-43 — Dectin-1 recognition of β-glucan; CARD9 → NF-κB. Source for: Dectin-1 rules.

## STAT signaling

- **Murray (2007)** *Curr Opin Pharmacol* 7(4):425-431 — STAT1 (IFN-γ, M1), STAT6 (IL-4, M2a), STAT3 (IL-10, M2c) selective activation. Source for: TF selectivity per cytokine in `transcription/`.

## Metabolism

- **O'Neill & Pearce (2016)** *J Exp Med* 213(1):15-23 — M1 = glycolysis (Warburg-like, succinate accumulation), M2 = OXPHOS / FAO. Source for: binary metabolic mode coupled to polarization in `metabolism/`.

## Phagocytosis & effector functions

- **Mosser & Edwards (2008)** *Nat Rev Immunol* (above) — classical effector functions: phagocytosis, ROS burst, NO, IL-12, TNF-α. Source for: effector rates in `effectors/`.
- **Underhill & Goodridge (2012)** *Nat Rev Immunol* 12(7):492-502 — phagocytosis mechanism: receptor engagement (FcγR for opsonized, Dectin-1 for fungal) → cup formation → engulfment. Source for: cup-state morphology + receptor-gated phagocytosis logic.
- **Forman & Torres (2002)** *Am J Respir Crit Care Med* 166(suppl):S4-S8 — ROS burst kinetics post-phagocytosis. Source for: ROS production rate.

## Cytokine secretion

- **Arango Duque & Descoteaux (2014)** *Front Immunol* 5:491 — quantitative cytokine secretion profiles (TNF-α, IL-6, IL-1β, IL-10, TGF-β) per polarization state. Source for: secretion rates in `effectors/cytokines.py`.

## Motility / chemotaxis

- **Lauffenburger & Linderman (1993)** *Receptors: Models for Binding, Trafficking, and Signaling* — classical biased-random-walk model for chemotaxis: drift = `χ * ∇log(L) / (1 + L/K)`. Source for: `motility/chemotaxis.py`.
- **Deshmane et al. (2009)** *J Interferon Cytokine Res* 29(6):313-326 — CCL2/CCR2 axis for monocyte/macrophage recruitment. Source for: CCL2 as the primary chemokine in v0.1.

## Morphology

- **Adams & Hamilton (1984)** *Annu Rev Immunol* 2:283-318 — classical descriptions of macrophage morphology: round (resting), spread/dendritic (activated), cup-shaped (phagocytosing). Source for: discrete morphology states.

## What is *not* yet grounded (v0.1 limits)

- Per-cell-cycle dynamics (proliferation, apoptosis) — out of scope.
- Trained immunity / epigenetic memory (Netea et al. 2016) — v0.3+.
- Macrophage subsets (Kupffer, microglia, alveolar) — v0.3+.
- Quantitative parameter calibration to specific cell lines — v0.2.

If you add a behavior, add the citation here and `# Ref:` it in code. A behavior without a citation is a fabrication and will be rejected by review.

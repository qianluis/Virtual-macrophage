"""Visualization helpers for the golden-standard experiments.

Each experiment's `run()` returns a dict of snapshots; these helpers turn them
into matplotlib figures saved to experiments/figures/. Optional dependency —
the deterministic core never imports matplotlib, so the package works headless
without it.

Run all plots:
    python -m macrophage.viz
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

FIG_DIR = Path(__file__).resolve().parents[2] / "experiments" / "figures"


def _ensure_fig_dir() -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    return FIG_DIR


def plot_polarization_timecourse(
    snapshots: list[tuple[int, dict[str, Any]]],
    title: str,
    filename: str,
) -> Path:
    import matplotlib.pyplot as plt

    times = [t for t, _ in snapshots]
    scores = [s["polarization_score"] for _, s in snapshots]
    labels = [s["polarization_label"] for _, s in snapshots]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(times, scores, "o-", linewidth=2, markersize=7, color="#2c7fb8")
    ax.axhline(0.3, color="#d62728", linestyle="--", alpha=0.5, label="M1 threshold")
    ax.axhline(-0.3, color="#1b9e77", linestyle="--", alpha=0.5, label="M2 threshold")
    ax.axhline(0, color="gray", linestyle=":", alpha=0.4)
    ax.set_xlabel("time (ticks, ~1 min each)")
    ax.set_ylabel("polarization score")
    ax.set_title(title)
    ax.set_ylim(-1.1, 1.1)
    for t, s, lab in zip(times, scores, labels):
        ax.annotate(lab, (t, s), textcoords="offset points", xytext=(0, 10),
                    ha="center", fontsize=8)
    ax.legend(loc="lower right")
    fig.tight_layout()
    out = _ensure_fig_dir() / filename
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_tf_panel(
    snapshots: list[tuple[int, dict[str, Any]]],
    title: str,
    filename: str,
) -> Path:
    import matplotlib.pyplot as plt

    times = [t for t, _ in snapshots]
    tfs = ["nfkb", "irf3", "stat1", "stat6", "stat3"]
    colors = ["#d62728", "#ff7f0e", "#9467bd", "#1b9e77", "#17becf"]

    fig, ax = plt.subplots(figsize=(7, 4))
    for tf, c in zip(tfs, colors):
        vals = [s["tfs"][tf] for _, s in snapshots]
        ax.plot(times, vals, "o-", label=tf.upper(), color=c, linewidth=2)
    ax.set_xlabel("time (ticks)")
    ax.set_ylabel("TF activity")
    ax.set_title(title)
    ax.set_ylim(-0.05, 1.1)
    ax.legend(ncol=5, fontsize=8)
    fig.tight_layout()
    out = _ensure_fig_dir() / filename
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_cytokine_output(
    snapshots: list[tuple[int, dict[str, Any]]],
    title: str,
    filename: str,
) -> Path:
    import matplotlib.pyplot as plt

    times = [t for t, _ in snapshots]
    cyto_keys = ["tnf_alpha", "il6", "il1b", "il10", "tgfb"]
    colors = ["#d62728", "#ff7f0e", "#9467bd", "#1b9e77", "#17becf"]

    fig, ax = plt.subplots(figsize=(7, 4))
    for k, c in zip(cyto_keys, colors):
        vals = [s["cytokines_secreted_total"][k] for _, s in snapshots]
        ax.plot(times, vals, "o-", label=k, color=c, linewidth=2)
    ax.set_xlabel("time (ticks)")
    ax.set_ylabel("cumulative secretion (ng/mL)")
    ax.set_title(title)
    ax.legend(ncol=5, fontsize=8)
    fig.tight_layout()
    out = _ensure_fig_dir() / filename
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_chemotaxis_trajectory(trajectories: list[np.ndarray], filename: str) -> Path:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 6))
    for traj in trajectories:
        ax.plot(traj[:, 0], traj[:, 1], "-", alpha=0.4, color="#2c7fb8", linewidth=0.8)
        ax.plot(*traj[0], "o", color="#2c7fb8", markersize=3)
        ax.plot(*traj[-1], "v", color="#d62728", markersize=5)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.set_title("Chemotaxis trajectories (CCL2 +x gradient)\nblue=start, red=end")
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    fig.tight_layout()
    out = _ensure_fig_dir() / filename
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main() -> None:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "experiments"))

    from experiments import lps_response, il4_response, chemotaxis  # noqa: E402

    print("plotting LPS response...")
    r = lps_response.run()
    plot_polarization_timecourse(r["snapshots"], "LPS (100 ng/mL) → M1 polarization",
                                 "lps_polarization.png")
    plot_tf_panel(r["snapshots"], "LPS: transcription factor activation",
                  "lps_tfs.png")
    plot_cytokine_output(r["snapshots"], "LPS: cytokine secretion (cumulative)",
                         "lps_cytokines.png")

    print("plotting IL-4 response...")
    r = il4_response.run()
    plot_polarization_timecourse(r["snapshots"], "IL-4 (20 ng/mL) → M2a polarization",
                                 "il4_polarization.png")
    plot_tf_panel(r["snapshots"], "IL-4: transcription factor activation\n(STAT3 = emergent autocrine IL-10 loop)",
                  "il4_tfs.png")
    plot_cytokine_output(r["snapshots"], "IL-4: cytokine secretion (cumulative)",
                         "il4_cytokines.png")

    print("plotting chemotaxis trajectories...")
    trajs = [chemotaxis.run(seed=s, n_ticks=200, with_gradient=True)["trajectory"]
             for s in range(12)]
    plot_chemotaxis_trajectory(trajs, "chemotaxis_trajectories.png")

    print(f"\nfigures written to {FIG_DIR}/")


if __name__ == "__main__":
    main()

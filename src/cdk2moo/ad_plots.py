"""Figures for the applicability-domain audit."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cdk2moo import config
from cdk2moo.ga_plots import INK, MUTED, GRID

SPLIT_COLOR = {"random": "#2a78d6", "scaffold": "#eb6834", "paper": "#1baf7a", "pooled": INK}


def _axes_style(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9)


def plot_reliability(stage4, path, title):
    """
    Mean absolute error of the forest against each similarity measure, by split,
    in eight equal-count bins of the measure, so every point rests on the same
    number of test predictions (rows; the same molecule recurs across seeds, so the
    number of distinct molecules is lower, see results/07_reliability_bins.csv).
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2), facecolor="#fcfcfb")
    panels = [("max_tanimoto_to_train", "raw ECFP4 max-Tanimoto to training set"),
              ("size_pctile", "size-conditioned similarity percentile")]
    for ax, (col, xlabel) in zip(axes, panels):
        _axes_style(ax)
        for name in ["random", "scaffold", "paper", "pooled"]:
            d = stage4 if name == "pooled" else stage4[stage4["split"] == name]
            edges = np.unique(np.quantile(d[col], np.linspace(0, 1, 9)))
            which = np.clip(np.digitize(d[col], edges[1:-1]), 0, len(edges) - 2)
            xs, ys = [], []
            for k in range(len(edges) - 1):
                g = d[which == k]
                xs.append(g[col].median())
                ys.append((g["pred_rf"] - g["y_true"]).abs().mean())
            ax.plot(xs, ys, marker="o", markersize=5, linewidth=2.4 if name == "pooled" else 1.6,
                    color=SPLIT_COLOR[name], label=name, zorder=3)
        ax.set_xlabel(xlabel, color=INK, fontsize=10)
        ax.set_ylabel("mean |predicted - measured| pActivity", color=INK, fontsize=10)
        ax.set_ylim(0, None)
    axes[0].legend(frameon=False, fontsize=9, title="Stage 4 test set", title_fontsize=9)
    fig.suptitle(title, color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_generated_on_domain(stage4, front, path, title):
    """
    Predicted pActivity against raw similarity. Grey points: Stage 4 test molecules
    (where we KNOW the error). Circles: the Stage 6 Pareto front, filled if it carries
    a PAINS or Brenk alert. The shaded box is where the front sits but Stage 4 has too
    few distinct molecules to estimate error at all.
    """
    fig, ax = plt.subplots(figsize=(10, 6.4), facecolor="#fcfcfb")
    _axes_style(ax)
    sample = stage4.drop_duplicates(["split", "row"])
    ax.scatter(sample["max_tanimoto_to_train"], sample["pred_rf"], s=6, color="#9a9a94",
               alpha=0.35, linewidth=0, label="Stage 4 out-of-sample test molecules")
    plain = front[~front["any_alert"]]
    flagged = front[front["any_alert"]]
    ax.scatter(plain["max_tanimoto"], plain["pred"], s=50, facecolor="none", edgecolor="#2a78d6",
               linewidth=1.5, label="Pareto front, no PAINS/Brenk alert")
    ax.scatter(flagged["max_tanimoto"], flagged["pred"], s=50, color="#eb6834", edgecolor="#fcfcfb",
               linewidth=0.8, label="Pareto front, with alert")
    for edge in config.AD_RAW_EDGES[1:-1]:
        ax.axvline(edge, color=MUTED, linewidth=0.8, linestyle=":")
    ax.axhline(8, color=MUTED, linewidth=0.8, linestyle=":")
    ax.fill_between([0, 0.6], 8, ax.get_ylim()[1] if False else 9.6, color="#eda100", alpha=0.10,
                    linewidth=0, label="high prediction, weak support: little or no Stage 4 data")
    ax.set_xlim(0.1, 1.0)
    ax.set_ylim(3.5, 9.6)
    ax.set_xlabel("raw ECFP4 (Morgan r=2, 2048-bit) max-Tanimoto to training set", color=INK, fontsize=10)
    ax.set_ylabel("predicted pActivity (real surrogate)", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle(title, color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=150)
    plt.close(fig)

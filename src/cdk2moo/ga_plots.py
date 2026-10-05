"""
Figures for the GA stage: trajectories per arm, and final-generation spread.

Colours are categorical slots 1-4 of the validated palette, always assigned to
the same arm. Two of the four sit below 3:1 contrast on white, so every arm
also has its own line style and a direct label at the line end.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ARM_STYLE = {
    "multi_real":      {"color": "#2a78d6", "ls": "-",  "label": "multi_real (activity x QED x SA)"},
    "multi_scrambled": {"color": "#eb6834", "ls": "--", "label": "multi_scrambled (scrambled surrogate)"},
    "activity_only":   {"color": "#1baf7a", "ls": ":",  "label": "activity_only"},
    "druglike_only":   {"color": "#eda100", "ls": "-.", "label": "druglike_only (QED x SA)"},
}
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e3e2dd"

PANELS = [
    ("best_pred_real", "best predicted pActivity (real surrogate)"),
    ("median_pred_real", "median predicted pActivity (real surrogate)"),
    ("median_max_tanimoto", "median max ECFP4 Tanimoto to training set"),
    ("median_size_pctile", "median size-conditioned similarity percentile"),
    ("median_qed", "median QED"),
    ("median_sa", "median SA score (lower = easier)"),
    ("median_heavy_atoms", "median heavy-atom count"),
    ("median_followed_pred", "median prediction of the surrogate the arm follows"),
]


def _style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def plot_trajectories(per_run, path, title):
    """
    Eight panels over generation (the last shows, for each arm, the surrogate it
    is actually optimising; druglike_only follows none, so has no line there). Each arm: line = median across GA seeds of the
    per-seed statistic, band = min to max across seeds (5 seeds: a range, not
    a confidence interval).
    """
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.5), facecolor="#fcfcfb")
    for ax, (col, ylabel) in zip(axes.ravel(), PANELS):
        _style_axes(ax)
        for arm, st in ARM_STYLE.items():
            sub = per_run[per_run["arm"] == arm]
            table = sub.pivot(index="generation", columns="seed", values=col)
            gens = table.index.to_numpy()
            ax.fill_between(gens, table.min(axis=1), table.max(axis=1),
                            color=st["color"], alpha=0.15, linewidth=0)
            ax.plot(gens, table.median(axis=1), color=st["color"], ls=st["ls"],
                    linewidth=2, label=st["label"])
        ax.set_xlabel("generation", color=MUTED, fontsize=9)
        ax.set_ylabel(ylabel, color=INK, fontsize=9)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle(title, color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_final_generation(final, path, title):
    """
    One dot per GA seed (that seed's final-generation population median), one
    column per arm, with a bar at the median over seeds. `final` has one row per
    arm x seed.
    """
    cols = [("median_pred_real", "predicted pActivity (real)"),
            ("median_max_tanimoto", "max ECFP4 Tanimoto to training"),
            ("median_size_pctile", "size-conditioned percentile"),
            ("median_qed", "QED"),
            ("median_sa", "SA score"),
            ("median_heavy_atoms", "heavy atoms")]
    fig, axes = plt.subplots(1, 6, figsize=(20, 4.5), facecolor="#fcfcfb")
    arms = list(ARM_STYLE)
    for ax, (col, ylabel) in zip(axes, cols):
        _style_axes(ax)
        for k, arm in enumerate(arms):
            vals = final.loc[final["arm"] == arm, col].to_numpy()
            jitter = np.linspace(-0.12, 0.12, len(vals))
            ax.scatter(k + jitter, vals, s=40, color=ARM_STYLE[arm]["color"],
                       edgecolor="#fcfcfb", linewidth=1.5, zorder=3)
            ax.hlines(np.median(vals), k - 0.25, k + 0.25, color=INK, linewidth=2, zorder=4)
        ax.set_xticks(range(len(arms)))
        ax.set_xticklabels([a.replace("_", "\n") for a in arms], fontsize=8, color=MUTED)
        ax.set_ylabel(ylabel, color=INK, fontsize=9)
    fig.suptitle(title, color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)

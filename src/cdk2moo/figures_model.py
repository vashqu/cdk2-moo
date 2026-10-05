"""Final figures 1, 2 and 4: the surrogate, the optimisation trajectories, and the reward-hacking view."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import textwrap

import numpy as np

from cdk2moo.ga_plots import ARM_STYLE, INK, MUTED, GRID
from cdk2moo.ad_plots import SPLIT_COLOR

LIMS = (3.4, 9.9)


def _style(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9)


def fig_surrogate_honesty(pred, metrics, seed, path, tag=""):
    """
    Figure 1. Predicted against measured pActivity for test molecules (one seed), random vs
    scaffold vs paper split in the top row, and the same forest trained on SCRAMBLED labels in
    the bottom row. Titles give R2 as median [25th, 75th percentile] over all ten split seeds.
    """
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), facecolor="#fcfcfb", sharex=True, sharey=True)
    for col, split in enumerate(["random", "scaffold", "paper"]):
        d = pred[(pred["seed"] == seed) & (pred["split"] == split)]
        for row, (column, model, color, name) in enumerate(
                [("pred_rf", "rf", SPLIT_COLOR[split], "random forest"),
                 ("pred_scrambled", "scrambled", "#9a9a94", "label-scrambled control")]):
            ax = axes[row, col]
            _style(ax)
            ax.scatter(d["y_true"], d[column], s=9, color=color, alpha=0.55, linewidth=0)
            ax.plot(LIMS, LIMS, color=INK, linewidth=1, linestyle="--")
            r2 = metrics[(metrics["split"] == split) & (metrics["model"] == model)]["r2"]
            q = np.percentile(r2, [25, 50, 75])
            ax.set_title(f"{split} split, {name}\nR2 {q[1]:.2f} [{q[0]:.2f}, {q[2]:.2f}] over 10 seeds",
                         color=INK, fontsize=10)
            ax.set_xlim(LIMS)
            ax.set_ylim(LIMS)
            if col == 0:
                ax.set_ylabel("predicted pActivity", color=INK)
            if row == 1:
                ax.set_xlabel("measured pActivity", color=INK)
    fig.suptitle(textwrap.fill(f"Predicted vs measured pActivity for held-out molecules (points: split seed {seed}; dashed: perfect prediction). {tag}", 125),
                 color=INK, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_trajectories(per_run, path, tag=""):
    """
    Figure 2. Four statistics against generation, mean and one standard deviation over the five GA
    seeds, for the four arms. Generation 0 = the starting training molecules; later generations
    describe molecules born after generation 0.
    """
    panels = [("median_pred_real", "median predicted pActivity (real surrogate)"),
              ("median_max_tanimoto", "median raw max ECFP4 Tanimoto to training"),
              ("median_qed", "median QED"), ("median_sa", "median SA score (lower = easier)")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), facecolor="#fcfcfb")
    for ax, (column, label) in zip(axes.ravel(), panels):
        _style(ax)
        for arm, st in ARM_STYLE.items():
            table = per_run[per_run["arm"] == arm].pivot(index="generation", columns="seed", values=column)
            mean, sd = table.mean(axis=1), table.std(axis=1)
            ax.fill_between(table.index, mean - sd, mean + sd, color=st["color"], alpha=0.18, linewidth=0)
            ax.plot(table.index, mean, color=st["color"], ls=st["ls"], linewidth=2, label=st["label"])
        ax.set_xlabel("generation", color=INK)
        ax.set_ylabel(label, color=INK, fontsize=9)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="center right")
    fig.suptitle(textwrap.fill(f"Median over generated molecules per generation (generation 0 = starting molecules); line = mean, band = +/- 1 SD over {per_run['seed'].nunique()} GA seeds. {tag}", 125),
                 color=INK, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_reward_hacking(gen, stage4, path, tag=""):
    """
    Figure 4. Each point is a distinct generated molecule: predicted pActivity against raw max
    ECFP4 Tanimoto to the training set, coloured by the generation it was born in (darker =
    later). The inset is the surrogate's measured error on real held-out molecules against the same
    similarity (Stage 4, pooled over splits and seeds, eight equal-count bins): where similarity is
    low, a high prediction is unvalidated.
    """
    edges = np.unique(np.quantile(stage4["max_tanimoto_to_train"], np.linspace(0, 1, 9)))
    which = np.clip(np.digitize(stage4["max_tanimoto_to_train"], edges[1:-1]), 0, len(edges) - 2)
    inset_x, inset_y = [], []
    for k in range(len(edges) - 1):
        g = stage4[which == k]
        inset_x.append(g["max_tanimoto_to_train"].median())
        inset_y.append((g["pred_rf"] - g["y_true"]).abs().mean())

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor="#fcfcfb", sharex=True, sharey=True)
    for ax, (arm, st) in zip(axes.ravel(), ARM_STYLE.items()):
        _style(ax)
        d = gen[gen["arm"] == arm].drop_duplicates("smiles").sort_values("birth_generation")
        sc = ax.scatter(d["max_tanimoto"], d["pred_real"], c=d["birth_generation"], cmap="Blues",
                        s=7, alpha=0.8, linewidth=0, vmin=0, vmax=50)
        ax.axvline(0.6, color=MUTED, linewidth=0.8, linestyle=":")
        ax.axhline(8, color=MUTED, linewidth=0.8, linestyle=":")
        ax.set_title(st["label"], color=INK, fontsize=10)
        ax.set_ylim(2.4, 9.0)             # leaves an empty strip under the data for the inset
        inset = ax.inset_axes([0.60, 0.08, 0.37, 0.22])
        inset.plot(inset_x, inset_y, color=INK, marker="o", markersize=3, linewidth=1.3)
        inset.set_ylim(0, 1.1)
        inset.set_title("real-molecule error vs similarity", fontsize=7, color=INK)
        inset.set_ylabel("mean |error|", fontsize=7, color=MUTED)
        inset.tick_params(labelsize=7, colors=MUTED)
        for side in ("top", "right"):
            inset.spines[side].set_visible(False)
    for ax in axes[1]:
        ax.set_xlabel("raw max ECFP4 Tanimoto to training set", color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel("predicted pActivity (real surrogate)", color=INK)
    fig.colorbar(sc, ax=axes, label="generation born", shrink=0.6, pad=0.02)
    fig.suptitle(textwrap.fill(f"Predicted activity against raw ECFP4 similarity to the training set, one point per distinct generated molecule (dotted: similarity 0.6, prediction 8). {tag}", 125),
                 color=INK, x=0.02, ha="left", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)

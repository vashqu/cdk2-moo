"""Figures for the Pareto stage: three 2-D views of the 3-objective trade-off, and a structure grid."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw

from cdk2moo.ga_plots import ARM_STYLE, INK, MUTED, GRID

VIEWS = [("qed", "pred", "QED (higher = better)", "predicted pActivity, real surrogate"),
         ("sa", "pred", "SA score (lower = better)", "predicted pActivity, real surrogate"),
         ("qed", "sa", "QED (higher = better)", "SA score (lower = better)")]


def plot_fronts(by_arm, front, actives, reference, path, title):
    """
    by_arm    {arm: DataFrame of unique generated molecules, columns pred, qed, sa}
    front     DataFrame of the pooled Pareto front (same columns)
    actives   known actives: measured pActivity in `pred`
    reference random training molecules: measured pActivity in `pred`
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.6), facecolor="#fcfcfb")
    for ax, (x, y, xlabel, ylabel) in zip(axes, VIEWS):
        ax.set_facecolor("#fcfcfb")
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=9)

        ax.scatter(reference[x], reference[y], s=14, marker="s", color="#9a9a94",
                   alpha=0.5, linewidth=0, label="random training molecules (measured)")
        ax.scatter(actives[x], actives[y], s=26, marker="D", facecolor="none",
                   edgecolor=INK, linewidth=0.9, label="known actives, top decile (measured)")
        for arm, st in ARM_STYLE.items():
            d = by_arm[arm]
            ax.scatter(d[x], d[y], s=5, color=st["color"], alpha=0.3, linewidth=0,
                       label=st["label"])
        ax.scatter(front[x], front[y], s=46, facecolor="none", edgecolor=INK,
                   linewidth=1.3, label="pooled Pareto front (generated)")
        ax.set_xlabel(xlabel, color=INK, fontsize=10)
        ax.set_ylabel(ylabel, color=INK, fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=9,
               markerscale=2.2)
    fig.suptitle(title, color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.11, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def draw_structures(smiles, legends, path):
    """Render molecules as a grid image with a caption under each."""
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(520, 380), legends=legends)
    img.save(str(path))

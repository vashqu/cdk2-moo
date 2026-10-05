"""Final figure 5: the orthogonal check, docking scores, size, pose validity and hinge contacts."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cdk2moo.ga_plots import ARM_STYLE, INK, MUTED, GRID

SETS = ["known_actives", "random_training", "decoys", "multi_real", "multi_scrambled",
        "activity_only", "druglike_only", "pareto_front"]
COLORS = {"known_actives": INK, "random_training": "#9a9a94", "decoys": "#52514e",
          "multi_real": ARM_STYLE["multi_real"]["color"], "multi_scrambled": ARM_STYLE["multi_scrambled"]["color"],
          "activity_only": ARM_STYLE["activity_only"]["color"], "druglike_only": ARM_STYLE["druglike_only"]["color"],
          "pareto_front": "#e87ba4"}


def _style(ax):
    ax.set_facecolor("#fcfcfb")
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(colors=MUTED, labelsize=9)


def fig_orthogonal(df, path, seed, pb_note=""):
    """
    Figure 5. df has one row per docked molecule: set, vina_score, n_heavy_atoms, pb_pass, hinge_contact.
    (a) Vina score per set, bar = median. (b) score against size, with the line fitted on random + decoy
    molecules. (c) share of poses passing all the PoseBusters checks that were audited (`pb_note` says which audit). (d) share of poses with a ligand N or O
    within 3.5 A of a hinge backbone atom (Glu81 O, Leu83 N or O): a distance-only proximity measure, not a verified hydrogen bond.
    """
    rng = np.random.default_rng(seed)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor="#fcfcfb")
    ax = axes[0, 0]
    _style(ax)
    for k, name in enumerate(SETS):
        d = df[df["set"] == name]
        ax.scatter(d["vina_score"], np.full(len(d), k) + rng.uniform(-0.28, 0.28, len(d)), s=9,
                   color=COLORS[name], alpha=0.6, linewidth=0)
        ax.vlines(d["vina_score"].median(), k - 0.4, k + 0.4, color=INK, linewidth=2)
    ax.set_yticks(range(len(SETS)))
    ax.set_yticklabels(SETS, fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel("Vina score, kcal/mol (more negative = better); bar = median", color=INK)
    ax.set_title("(a) known actives, decoys and random molecules overlap", color=INK, fontsize=10)

    ax = axes[0, 1]
    _style(ax)
    for name in SETS:
        d = df[df["set"] == name]
        ax.scatter(d["n_heavy_atoms"], d["vina_score"], s=12, color=COLORS[name], alpha=0.6, linewidth=0, label=name)
    ref = df[df["set"].isin(["random_training", "decoys"])]
    slope, intercept = np.polyfit(ref["n_heavy_atoms"], ref["vina_score"], 1)
    xs = np.array([df["n_heavy_atoms"].min(), df["n_heavy_atoms"].max()])
    ax.plot(xs, intercept + slope * xs, color=INK, linestyle="--", linewidth=1.4, label="fit on random + decoys")
    ax.set_xlabel("heavy atoms", color=INK)
    ax.set_ylabel("Vina score, kcal/mol", color=INK)
    ax.set_title("(b) score against size", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper right")

    for ax, column, title in [(axes[1, 0], "pb_pass", f"(c) docked poses passing all audited PoseBusters checks (%){pb_note}"),
                              (axes[1, 1], "hinge_contact", "(d) poses with a ligand N/O within 3.5 A of the hinge backbone (%), distance only")]:
        _style(ax)
        rates = [100 * df[df["set"] == s][column].mean() for s in SETS]
        ax.barh(range(len(SETS)), rates, color=[COLORS[s] for s in SETS], height=0.6)
        ax.set_yticks(range(len(SETS)))
        ax.set_yticklabels(SETS, fontsize=9, color=INK)
        ax.invert_yaxis()
        ax.set_xlim(0, 105)
        for k, r in enumerate(rates):
            ax.text(r + 1, k, f"{r:.0f}%", va="center", fontsize=9, color=INK)
        ax.set_title(title, color=INK, fontsize=10)
    fig.suptitle("Figure 5. Docking as an orthogonal check (4KD1, Vina seed 42): scores do not separate actives; poses are valid; "
                 "activity-driven molecules reach the hinge", color=INK, x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_docking_group(g, group, role_text, campaign_id, scope, path):
    """
    Vina scores for one experiment group, from completed jobs only. Left: the score of every distinct molecule per set (bar = median; the number of
    completed molecules is in each label). Right: score against heavy-atom count, with the straight-line fit over ALL molecules of the group, because
    Vina scores improve with size. The comparator actives and their role (held out or in-sample) are stated in the title.
    """
    order = [s for s in ["training_top_decile", "heldout_actives", "random_training", "random_remaining", "decoys", "multi_real", "multi_scrambled",
                         "activity_only", "druglike_only", "pareto_front"] if f"{group}:{s}" in set(g["set"])]
    palette = {"training_top_decile": INK, "heldout_actives": INK, "random_training": "#9a9a94", "random_remaining": "#9a9a94", "decoys": "#52514e",
               "multi_real": ARM_STYLE["multi_real"]["color"], "multi_scrambled": ARM_STYLE["multi_scrambled"]["color"],
               "activity_only": ARM_STYLE["activity_only"]["color"], "druglike_only": ARM_STYLE["druglike_only"]["color"], "pareto_front": "#e87ba4"}
    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.6), facecolor="#fcfcfb")
    ax = axes[0]
    _style(ax)
    for k, name in enumerate(order):
        d = g[g["set"] == f"{group}:{name}"]
        ax.scatter(d["vina_score"], np.full(len(d), k) + rng.uniform(-0.28, 0.28, len(d)), s=9, color=palette[name], alpha=0.6, linewidth=0)
        ax.vlines(d["vina_score"].median(), k - 0.4, k + 0.4, color=INK, linewidth=2)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{n} (n={int((g['set'] == f'{group}:{n}').sum())})" for n in order], fontsize=9, color=INK)
    ax.invert_yaxis()
    ax.set_xlabel("Vina score, kcal/mol (more negative = better); bar = median; points = completed molecules", color=INK, fontsize=9)
    ax = axes[1]
    _style(ax)
    for name in order:
        d = g[g["set"] == f"{group}:{name}"]
        ax.scatter(d["n_heavy_atoms"], d["vina_score"], s=12, color=palette[name], alpha=0.6, linewidth=0, label=name)
    slope, intercept = np.polyfit(g["n_heavy_atoms"], g["vina_score"], 1)
    xs = np.array([g["n_heavy_atoms"].min(), g["n_heavy_atoms"].max()])
    ax.plot(xs, intercept + slope * xs, color=INK, linestyle="--", linewidth=1.4, label=f"linear fit, all molecules ({slope:+.3f}/atom)")
    ax.set_xlabel("heavy atoms", color=INK)
    ax.set_ylabel("Vina score, kcal/mol", color=INK)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper right")
    fig.suptitle(f"Docking into 4KD1 (Vina), group {group}: {role_text}.  Campaign {campaign_id}, {scope} policy, Vina seed 42; distinct molecules, no pooling across groups or policies",
                 color=INK, x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_pose_checks(pb_by_set, hinge_by_set, tag, path):
    """
    Left: for every docked set, the share of its molecules whose pose passed all 22 PoseBusters checks, failed at least one, was unevaluable, or has no audit
    result (denominator = molecules in the set). Right: the share of the set's molecules whose docked pose has a ligand N or O within 3.5 A of the hinge backbone
    (distance only; molecules without a pose count as no contact; the denominator is again the set).
    """
    fig, axes = plt.subplots(1, 2, figsize=(15, 7), facecolor="#fcfcfb", sharey=True)
    names = list(pb_by_set["set"])
    ax = axes[0]
    _style(ax)
    left = np.zeros(len(pb_by_set))
    for column, color in [("pass", "#2a78d6"), ("fail", "#e34948"), ("unevaluable", "#eda100"), ("no_audit_result", "#9a9a94")]:
        share = 100 * pb_by_set[column].to_numpy() / pb_by_set["n"].to_numpy()
        ax.barh(range(len(names)), share, left=left, color=color, height=0.6, label=column)
        left += share
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels([f"{n} (n={int(k)})" for n, k in zip(names, pb_by_set["n"])], fontsize=8, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of molecules in the set", color=INK)
    ax.set_title("PoseBusters, 22 required checks", color=INK, fontsize=10)
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.12))
    ax = axes[1]
    _style(ax)
    shares = 100 * hinge_by_set.set_index("set").loc[names, "share_of_set"].to_numpy()
    ax.barh(range(len(names)), shares, color="#1baf7a", height=0.6)
    for k, v in enumerate(shares):
        ax.text(v + 1, k, f"{v:.0f}%", va="center", fontsize=8, color=INK)
    ax.set_xlim(0, 105)
    ax.set_xlabel("% of molecules in the set with a ligand N/O within 3.5 A of the hinge backbone (distance only)", color=INK, fontsize=9)
    ax.set_title("Hinge proximity", color=INK, fontsize=10)
    fig.suptitle(f"Docked-pose checks. {tag}", color=INK, x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)

#!/usr/bin/env python3
"""
Stage 11a: the campaign's final figures for one policy scope, collected in figures/final/.

Reads : results/04_*, 05_ga_per_run.csv, 05_ga_populations_scored.csv, 06_pareto_*.png (figures/), 07_stage4_ad.csv, 08e_*, 08f/08g by-set tables
Writes: figures/final/fig1_surrogate.png            predicted vs measured pActivity, three splits, with the label-scrambled control
        figures/final/fig2_trajectories.png          four GA arms over 50 generations (mean +/- SD over 5 seeds)
        figures/final/fig3_pareto_fronts.png, fig3b_pareto_structures.png   copied from stage 6
        figures/final/fig4_similarity_vs_prediction.png   predicted activity vs similarity, coloured by generation, with the surrogate's measured error inset
        figures/final/fig5_docking_<group>.png      Vina scores per set for groups P, M, S (copied from stage 8e)
        figures/final/fig6_pose_checks.png          PoseBusters and hinge proximity per set, with denominators

Every figure is drawn from saved campaign results; titles state the quantity, the policy scope and the sample sizes, not a conclusion. Existing figures are
never replaced (pass --overwrite to redo).

Run (inside a policy scope, after the analysis stages):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/11_make_figures.py
"""

import argparse
import shutil

import pandas as pd

from cdk2moo import config
from cdk2moo.figures_model import fig_surrogate_honesty, fig_trajectories, fig_reward_hacking
from cdk2moo.figures_validation import fig_pose_checks


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()
    R, F = config.RESULTS_DIR, config.FIGURES_DIR
    out = F / "final"
    out.mkdir(exist_ok=True)
    tag = f"Campaign {config.CAMPAIGN}, {config.SCOPE} policy scope."
    targets = {n: out / n for n in ["fig1_surrogate.png", "fig2_trajectories.png", "fig3_pareto_fronts.png", "fig3b_pareto_structures.png",
                                    "fig4_similarity_vs_prediction.png", "fig6_pose_checks.png"]}
    for group in "PMS":
        targets[f"fig5_docking_{group}.png"] = out / f"fig5_docking_{group}.png"
    clash = [str(p) for p in targets.values() if p.exists()]
    if clash and not args.overwrite:
        raise SystemExit(f"figures exist and are not replaced: {clash[:3]}...")

    pred = pd.read_csv(R / "04_test_predictions.csv")
    metrics = pd.read_csv(R / "04_surrogate_metrics.csv")
    fig_surrogate_honesty(pred, metrics, config.RANDOM_SEED, targets["fig1_surrogate.png"], tag + " Surrogate trained on each split's training set; 10 split seeds in the titles (median [quartiles]).")
    fig_trajectories(pd.read_csv(R / "05_ga_per_run.csv"), targets["fig2_trajectories.png"], tag + " Primary experiment (surrogate trained on all 2,016 molecules).")
    shutil.copy(F / "06_pareto_fronts.png", targets["fig3_pareto_fronts.png"])
    shutil.copy(F / "06_pareto_structures.png", targets["fig3b_pareto_structures.png"])
    pop = pd.read_csv(R / "05_ga_populations_scored.csv", usecols=["smiles", "arm", "pred_real", "max_tanimoto", "birth_generation"])
    fig_reward_hacking(pop[pop["birth_generation"] > 0], pd.read_csv(R / "07_stage4_ad.csv"), targets["fig4_similarity_vs_prediction.png"],
                       tag + " Primary experiment; inset: Stage 4 out-of-sample error.")
    for group in "PMS":
        shutil.copy(F / f"08_docking_{group}.png", targets[f"fig5_docking_{group}.png"])
    pb = pd.read_csv(R / "posebusters_audit" / "posebusters_by_set.csv")
    hinge = pd.read_csv(R / "08g_hinge_by_set.csv")
    fig_pose_checks(pb, hinge, tag + " Vina seed 42 poses.", targets["fig6_pose_checks.png"])
    for p in sorted(targets.values()):
        print(f"  {p.relative_to(config.CAMPAIGN_ROOT)}")


if __name__ == "__main__":
    main()

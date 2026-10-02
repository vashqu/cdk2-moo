# CDK2 multi-objective molecular optimization

**Question.** When a molecular optimizer is driven by a machine-learned activity surrogate under
drug-likeness constraints, does improvement in the objective reflect progress toward the target, or
exploitation of the surrogate's extrapolation error?

A negative or null result is a valid outcome here. Every judgment call is logged in
[`decisions.md`](decisions.md) (D-01 to D-32); project rules and known confounds are in
[`CLAUDE.md`](CLAUDE.md).

## Pipeline

Run in order from the repo root, in the `cdk2` conda env. Stages talk only through files.

| Script | What it does | Main outputs |
|---|---|---|
| `01_fetch_chembl.py` | Download CDK2 (CHEMBL301) activities, cached | `data/raw/` |
| `02_curate.py` | Filter, standardize, deduplicate by structure: 2,016 molecules | `data/processed/cdk2_ic50_curated.csv` |
| `03_features_and_splits.py` | ECFP4 fingerprints; random, scaffold and paper splits, 10 seeds | `ecfp4.npy`, `cdk2_splits.csv` |
| `04_train_surrogate.py` | Random forest, 3 splits x 10 seeds, with scrambled-label, mean and size-only controls | `results/04_*` |
| `05_run_ga.py` | Genetic algorithm, 4 arms x 5 seeds, surrogate trained on all 2,016 | `results/05_ga_populations.csv` |
| `05b_size_conditioned_similarity.py` | Add a size-conditioned similarity percentile | `05_ga_populations_scored.csv` |
| `05c_summarize_ga.py` | Trajectories, final-generation tables, paired contrasts | `figures/05_ga_*.png` |
| `05d_ga_descriptive_checks.py` | Diversity collapse and size-boundary piling | `results/05_descriptive_checks.txt` |
| `05e_positive_control_setup.py` | Define held-out-actives controls (molecule and scaffold level); test the surrogate on them | `results/05e_*` |
| `05f_positive_control_ga.py` | Same GA with the held-out molecules removed from training (also runs the scaffold-trained robustness GA) | `results/05f_*_populations.csv` |
| `05g_positive_control_eval.py` | Did the GA find chemistry near the held-out actives? | `figures/05g_positive_control.png` |
| `05h_scaffold_surrogate_robustness.py` | Do the Stage 5 claims hold for a scaffold-trained surrogate? | `results/05h_robustness_summary.csv` |
| `08_prepare_receptor.py` | Prepare the 4KD1 receptor for Vina (protein only, hydrogens, PDBQT, search box) | `data/structures/receptor*` |
| `08b_redock_native_ligand.py` | Redock dinaciclib, RMSD to the crystal pose (0.65 A, 3 seeds) | `results/08b_redocking.json` |
| `08c_select_docking_sets.py` | Choose 8 sets to dock (seeded): known actives, random, decoys, four GA arms, Pareto front | `results/08c_docking_sets.csv` |
| `08d_dock_sets.py` | Dock 788 molecules with Vina (resumable, ~105 min); 30 replicated over 3 seeds | `results/08d_docking_scores.csv` |
| `08e_docking_analysis.py` | Do scores separate actives from decoys? Size control, H3 contrast, noise | `figures/08_docking_scores.png` |
| `08f_posebusters.py` | PoseBusters validity of the docked poses | `figures/08_posebusters.png` |
| `08g_hinge_contacts.py` | Do poses make the hinge hydrogen bonds of an ATP-site inhibitor? | `results/08g_hinge_contacts.csv` |
| `10_run_constrained_ga.py` | H5 mitigation arm: the GA with a hard similarity floor (0.4 to 0.7), 120 runs | `results/10_constrained_*_populations.csv` |
| `10b_analyze_tradeoff.py` | What the floor costs in predicted activity and buys in recovery of held-out actives | `figures/10_tradeoff.png` |
| `11_make_figures.py` | Final figures 1, 2, 4, 5 drawn from saved results; figure 3 copied from stage 6 | `figures/final/` |
| `11b_table1.py` | Table 1: per-set summary of the docked molecules | `results/table1.csv`, `figures/final/table1.md` |
| `06_pareto.py` | Pareto front over predicted pActivity, QED, SA | `figures/06_pareto_*.png` |
| `07_ad_reliability.py` | Surrogate error vs similarity to its training set | `figures/07_reliability_curves.png` |
| `07b_ad_generated_overlay.py` | Where the generated molecules sit on that curve | `results/07_generated_regions.csv` |
| `07c_structural_alerts.py` | PAINS / Brenk alert audit (recorded, not filtered) | `results/07_*_alerts.csv` |

Not done yet: the near-training-set mitigation arm (H5),
the random-ChEMBL reference sample, final figures and slides.

## Figures

All in `figures/final/`; each is drawn from saved results, nothing is recomputed.

| File | Shows |
|---|---|
| `fig1_surrogate_honesty.png` | Predicted vs measured pActivity on random, scaffold and paper splits, with the label-scrambled control |
| `fig2_trajectories.png` | Four GA arms over 50 generations: predicted activity, similarity to training, QED, SA (mean +/- SD, 5 seeds) |
| `fig3_pareto_fronts.png`, `fig3b_pareto_structures.png` | Pareto front over predicted pActivity, QED, SA with overlays; four front molecules drawn out |
| `fig4_reward_hacking.png` | Predicted activity vs similarity to training, coloured by generation, with the surrogate's measured error vs similarity as an inset |
| `fig5_orthogonal_validation.png` | Docking: scores per set, score vs size, PoseBusters pass rate, hinge contacts |
| `table1.md` | Per-set summary (n, activity, Vina, QED, SA, MW, cLogP, diversity, % PAINS, % PoseBusters pass, similarity) |

The random reference in figure 3 and Table 1 is a sample of training molecules, a stand-in for a random-ChEMBL
sample that has not been downloaded.

## Results so far (descriptive)

- **Surrogate.** Test R2 0.65 / 0.49 / 0.30 on random / scaffold / paper splits (median over 10 seeds);
  the scrambled-label and size-only controls score about zero.
- **Similarity alone is not exploitation.** In the GA, falling similarity to the training set happens
  just as much with a drug-likeness-only objective that never sees the surrogate.
- **Activity pressure finds known potent chemistry rather than extrapolating.** Activity-driven arms
  gain about +1.9 predicted pActivity over the drug-likeness baseline in every seed, but end no
  further from the training set (stable across five scaffold-split surrogates). H2 (as worded) is
  not supported. The multi-objective arm collapses onto few scaffolds in every split; the
  activity-only arm's collapse is fragile across splits (D-28).
- **Scrambled surrogate (H4).** The GA raises its own random surrogate's prediction (about +0.5), but
  the real surrogate's prediction for those molecules sits at the drug-likeness baseline (pooled
  difference -0.15, against +1.55 for the real-surrogate arm).
- **Reliability.** Surrogate error rises as similarity falls (mean absolute error 0.39 at raw similarity
  of 0.8 or above, 0.91 below 0.4). Stage 4 has almost no validation data where the prediction is 8 or more
  and similarity is below 0.6, so no reliability is claimed there.
- **Docking (H3).** The redocking control passes (0.65 A), but Vina does not separate known actives from
  property-matched decoys (AUC 0.53) or track measured activity (Spearman +0.02), so it cannot serve as an
  orthogonal check here and H3 is inconclusive. Activity-only molecules score better than known actives,
  the multi-objective molecules and Pareto front do not. Poses show the activity-driven search finds
  hinge-binding chemistry (96% of poses touch the hinge, against 44% for the drug-likeness-only arm).
- **Mitigation (H5).** Constraining the GA to stay near the training set does not lower predicted activity
  (it is flat or slightly higher), so H5's "attenuation" is not supported; the floor costs drug-likeness
  instead (multi-objective QED 0.91 to 0.72 at the strictest floor). Recovery of held-out actives improves
  only at the molecule level for the multi-objective arm, which keeps molecules near the actives' siblings;
  at the scaffold level there is no reliable gain.
- **Held-out actives.** The activity-only search recovers held-out potent molecules (sibling-supported)
  and partly recovers held-out scaffold families; the composite objective does not.

## Limitations

See the list in `CLAUDE.md`. In short: assay heterogeneity, no measurable noise floor, ECFP4 cannot
represent activity cliffs, GA edits are chemically naive, QED rewards resemblance to existing drugs, and
nothing is validated experimentally.

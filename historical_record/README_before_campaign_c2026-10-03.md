# CDK2 multi-objective molecular optimization

**Question.** When a molecular optimizer is driven by a machine-learned activity surrogate under
drug-likeness constraints, does improvement in the objective reflect progress toward the target, or
exploitation of the surrogate's extrapolation error?

A negative or null result is a valid outcome here. Every judgment call is logged in
[`decisions.md`](decisions.md) (append-only; later entries correct earlier ones rather than rewrite them);
project rules and known confounds are in [`CLAUDE.md`](CLAUDE.md); environment information is in
[`reproducibility/`](reproducibility/).

## Pipeline (run in this order, from the repo root, in the `cdk2` conda environment)

Stages communicate only through files. "Historical" outputs were produced before the repair pass described below
and are kept unchanged; scripts that would have replaced them now write to new, separately named locations.

| Script | What it does | Main outputs |
|---|---|---|
| `01_fetch_chembl.py` | Download CDK2 (CHEMBL301) activities, cached | `data/raw/` |
| `02_curate.py` | Filter, standardize, deduplicate by structure: 2,016 molecules | `data/processed/cdk2_ic50_curated.csv` |
| `03_features_and_splits.py` | ECFP4 fingerprints; random, scaffold and paper splits, 10 seeds | `ecfp4.npy`, `cdk2_splits.csv` |
| `04_train_surrogate.py` | Random forest, 3 splits x 10 seeds, with scrambled-label, mean and size-only controls | `results/04_*` |
| `05_run_ga.py` | Genetic algorithm, 4 arms x 5 seeds, surrogate trained on all 2,016. Output name depends on the settings; the primary file is protected (`--policy`, `--relaxed`, `--out`, `--overwrite`) | `results/05_ga_populations.csv` (+ `.manifest.json` for new runs) |
| `05b_size_conditioned_similarity.py` | Add a size-conditioned similarity percentile | `05_ga_populations_scored.csv` |
| `05c_summarize_ga.py` | Trajectories, final-generation tables, paired contrasts | `figures/05_ga_*.png` |
| `05d_ga_descriptive_checks.py` | Diversity collapse and size-boundary piling | `results/05_descriptive_checks.txt` |
| `05e_positive_control_setup.py` | Define held-out-actives controls (molecule and scaffold level); test the surrogate on them | `results/05e_*` |
| `05f_positive_control_ga.py` | Same GA with the held-out molecules removed from training (also the scaffold-trained robustness runs) | `results/05f_*_populations.csv` |
| `05g_positive_control_eval.py` | Proximity of generated molecules to the held-out actives (corrected definitions, `_v2` outputs) | `results/05g_positive_control_summary_v2.csv` |
| `05h_scaffold_surrogate_robustness.py` | Do the Stage 5 claims hold for scaffold-trained surrogates (5 splits)? | `results/05h_robustness_summary.csv` |
| `06_pareto.py` | Pareto front over predicted pActivity, QED, SA | `figures/06_pareto_*.png` |
| `07_ad_reliability.py` | Surrogate error vs similarity to its training set | `figures/07_reliability_curves.png` |
| `07b_ad_generated_overlay.py` | Where the generated molecules sit on that curve | `results/07_generated_regions.csv` |
| `07c_structural_alerts.py` | PAINS / Brenk alert audit (recorded, not filtered) | `results/07_*_alerts.csv` |
| `08_prepare_receptor.py` | Prepare the 4KD1 receptor for Vina | `data/structures/receptor*` |
| `08b_redock_native_ligand.py` | Redock dinaciclib, RMSD to the crystal pose (0.65 A, 3 seeds) | `results/08b_redocking.json` |
| `08c_select_docking_sets.py` | Choose 8 sets to dock (seeded) | `results/08c_docking_sets.csv` |
| `08d_dock_sets.py` | Dock 788 molecules with Vina. Resume-safe job store with a provenance manifest; writes to `--run-dir` | `results/dock_runs/seed42/` (new runs) |
| `08e_docking_analysis.py` | Do scores separate actives from decoys? Size control, noise | `figures/08_docking_scores.png` |
| `08f_posebusters.py` | PoseBusters audit of the docked poses, 22 required checks, three states | `results/posebusters_audit/` |
| `08g_hinge_contacts.py` | Ligand N/O within 3.5 A of the hinge backbone (distance only) | `results/08g_hinge_contacts.csv` |
| `10_run_constrained_ga.py` | H5 mitigation arm: GA with a hard similarity floor (0.4 to 0.7) | `results/10_constrained_*_populations.csv` |
| `10b_analyze_tradeoff.py` | Cost of the floor in predicted activity and its effect on proximity to held-out actives (`_v2` outputs) | `figures/10_tradeoff_v2.png` |
| `11_make_figures.py` | Final figures 1, 2, 4, 5; figure 3 copied from stage 6 | `figures/final/` |
| `11b_table1.py` | Table 1, per-set summary of the docked molecules | `results/table1*.csv`, `figures/final/table1*.md` |
| `12_audit_leakage_and_artifacts.py` | Read-only audits: split leakage, positive controls, docking coverage, PoseBusters coverage | `results/audit/` |

Tests: `python -m unittest discover tests` (standard library only).

Not done: the random-ChEMBL reference sample (the "random" sets are samples of the training molecules, a stand-in) and the
slide deck. The corrected GA variant (`--policy corrected`) exists but has **not** been run as a full campaign.

## Repair pass (2026-10-03): what changed, and what did not

Implementation repairs that preserve the experiment (historical results remain valid; behaviour was checked against saved
populations):
* `ga.py`: the similarity floor is now an explicit feasibility filter (infeasible molecules are removed before selection, the starting
  molecules must satisfy the floor, invalid thresholds raise). The legacy search reproduces the saved populations exactly (generations 0 to 3,
  seed 42, checked against `05_ga_populations.csv` and the tau = 0.6 constrained file).
* `05_run_ga.py` and the other GA entry points: the primary result filename is reserved; relaxed, smoke, and corrected runs get their own names;
  existing outputs are never replaced without `--overwrite`; each new run writes a manifest.
* `dock_jobs.py`, `docking.py`, `08d_dock_sets.py`: resume-safe docking (atomic per-job pose and record, hash manifest, grid-map validation,
  explicit retry and repair policy). The saved docking artifacts were audited and are complete (scores and poses agree exactly for every seed);
  they have no manifest and are treated as legacy.
* `posebusters_audit.py`, `08f_posebusters.py`: see the correction below.

Corrected analyses (new files; historical files kept):
* PoseBusters audit with all 22 required checks (`results/posebusters_audit/`), `fig5_orthogonal_validation_corrected_pb.png`, `table1_corrected_pb.*`.
  **774 of 787 poses (98.3%) pass all 22 checks, 13 fail, none is unevaluable.** The earlier "99% pass" came from a 21-check audit that dropped
  `internal_energy` (a check with missing values for some poses became object-typed and was silently excluded by a `dtype == bool` filter); six poses
  that audit passed fail that check.
* `05g_positive_control_summary_v2.csv`, `10_tradeoff_v2.*`: evaluation with explicit definitions. "Recovery" is now *proximity* (ECFP4 Tanimoto
  >= 0.6 to a held-out active); fingerprint identity, three exact-identity keys and distinct vs population-weighted counts are separate; the starting
  molecules no longer match themselves (their nearest-held-out share is 8% and 1%, not 0%); the old "chance rate" is only a reference-set prevalence;
  a neighbour's measured activity is an annotation of that neighbour, not a measurement of the generated molecule.
* `results/audit/`: leakage audits (below).

New experiment variants that need full runs before any conclusion (not run at scale):
* `--policy corrected`: every generated molecule is standardized like the training molecules (cleanup, largest fragment, neutralise, canonical
  tautomer), accepted only if closed-shell, deduplicated after standardization, and scored in that same form. The legacy search scores edited
  graphs as they are. The two are different experiments; only smoke tests of the mechanism were run.

Leakage and artifact audit (`results/audit/`): full standardized identity never crosses a split (duplicates were merged by InChIKey); stereoisomer-level
(connectivity) overlap does: median 23 test molecules in the random split, 3 in the paper split, 0 in the scaffold split; identical-fingerprint overlap is
27, 3 and 0. Four of the 204 molecule-level held-out actives have an identical-fingerprint twin left in training (none in the scaffold control). Every
held-out molecule in both controls shares at least one document with the remaining training molecules. `primary_document` is the first non-null document
met when records were aggregated, not a verified earliest publication, and the paper split is not a temporal split.

## Figures

All in `figures/final/`; drawn from saved results, nothing is recomputed.

| File | Shows |
|---|---|
| `fig1_surrogate_honesty.png` | Predicted vs measured pActivity on random, scaffold and paper splits, with the label-scrambled control |
| `fig2_trajectories.png` | Four GA arms over 50 generations: predicted activity, similarity to training, QED, SA (mean +/- SD, 5 seeds) |
| `fig3_pareto_fronts.png`, `fig3b_pareto_structures.png` | Pareto front over predicted pActivity, QED, SA with overlays; four front molecules drawn out |
| `fig4_reward_hacking.png` | Predicted activity vs similarity to training, coloured by generation, with the surrogate's measured error vs similarity as an inset |
| `fig5_orthogonal_validation.png` | Docking: scores, score vs size, PoseBusters pass rate (historical 21-check audit), distance-only hinge proximity. The corrected 22-check version is `fig5_orthogonal_validation_corrected_pb.png` |
| `table1.md`, `table1_corrected_pb.md` | Per-set summary (n, activity, Vina, QED, SA, MW, cLogP, diversity, % PAINS, % PoseBusters pass, similarity); the second uses the corrected audit |

The random reference in figure 3 and Table 1 is a sample of training molecules, a stand-in for a random-ChEMBL sample that has not been downloaded.
`pActivity` in Table 1 is measured for the known actives and random molecules and predicted for the rest.

## Results so far (descriptive; the repairs above do not change these conclusions)

- **Surrogate.** Test R2 0.65 / 0.49 / 0.30 on random / scaffold / paper splits (median over 10 seeds); the scrambled-label and size-only controls score about zero.
- **Similarity alone is not exploitation.** In the GA, falling similarity to the training set happens just as much with a drug-likeness-only objective that never sees the surrogate.
- **Activity pressure finds known potent chemistry rather than extrapolating.** Activity-driven arms gain about +1.9 predicted pActivity over the drug-likeness baseline in every
  seed but end no further from the training set (stable across five scaffold-split surrogates). H2 as worded is not supported. The multi-objective arm collapses onto few scaffolds
  in every split; the activity-only arm's collapse is fragile across splits (D-28).
- **Scrambled surrogate (H4).** The GA raises its own random surrogate's prediction (about +0.5), but the real surrogate's prediction for those molecules sits at the drug-likeness
  baseline (pooled difference -0.15, against +1.55 for the real-surrogate arm).
- **Reliability.** Surrogate error rises as similarity falls (mean absolute error 0.39 at raw similarity of 0.8 or above, 0.91 below 0.4). Stage 4 has almost no validation data where the
  prediction is 8 or more and similarity is below 0.6, so no reliability is claimed there.
- **Docking (H3).** The redocking control passes (0.65 A), but Vina does not separate known actives from property-matched decoys (AUC 0.53) or track measured activity (Spearman +0.02),
  so it cannot serve as an orthogonal check here and H3 is inconclusive. About 98% of docked poses pass the 22 PoseBusters checks, which is weak evidence because a docking engine
  places ligands without overlaps by construction. In 96% of activity-only poses a ligand N or O lies within 3.5 A of the hinge backbone (44% for the drug-likeness-only arm): a
  distance-only proximity measure, not a verified hydrogen bond.
- **Mitigation (H5).** A similarity floor does not lower predicted activity (it is flat or slightly higher), so H5's "attenuation" is not supported; it costs drug-likeness instead
  (multi-objective QED 0.91 to 0.72 at the strictest floor). The share of generated molecules *proximate* to held-out actives (ECFP4 >= 0.6) improves only at the molecule level for the
  multi-objective arm, which keeps molecules near the actives' siblings; at the scaffold level there is no reliable gain. This is fingerprint proximity, not recovered activity.
- **Held-out actives.** The activity-only search gets close (in fingerprint space) to held-out potent molecules (sibling-supported) and partly to held-out scaffold families; the composite
  objective does not. Exact standardized identity with a held-out active is 0.2% of activity-only molecules at the molecule level (about one molecule) and zero elsewhere.

## Things worth stating precisely

* **Elitism** keeps the best molecules *by scalar fitness*; it does not guarantee that every objective component improves monotonically (the composite can rise while a term falls).
* **"Newly proposed"** means not in the current population and not already proposed in the same batch. It does not mean globally unseen or novel: a proposal may have appeared in an
  earlier generation or run, or be a training molecule. No novelty filter is applied, since that would change the search.
* **The 20-39 heavy-atom window** is the central 5th-95th percentile of the training molecules' heavy-atom counts. The training set spans 5-78 heavy atoms, so about 10% of it lies outside.
* **Tautomer canonicalization is not always exhaustive.** For 23 of the 2,016 curated molecules (1.1%) RDKit's enumeration stopped at its transform limit, so two tautomers of one
  compound are not guaranteed to map to the same structure. RDKit's tautomer step also removes sp3 stereo at positions tautomerism could interconvert, so some stereoisomers (e.g.
  amino-acid enantiomers) standardize to the same string.
* **Vina** has an approximate rotatable-bond penalty but no rigorous treatment of binding entropy, no explicit water and a rigid receptor.

## Limitations

See the list in `CLAUDE.md`. In short: assay heterogeneity, no measurable noise floor, ECFP4 cannot tell apart molecules with identical fingerprints (e.g. stereoisomers) and
similarity-based models tend to smooth over activity cliffs, GA edits are chemically naive, QED rewards resemblance to existing drugs, and nothing is validated experimentally.
None of the repairs establishes biological activity, fixes the H3 docking problem, or validates H5 transfer; those remain open scientific questions.

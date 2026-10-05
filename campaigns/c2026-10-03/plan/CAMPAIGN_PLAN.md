# Campaign c2026-10-03: execution and analysis plan (frozen before any campaign stage ran)

This file is hashed and its hash recorded in `manifests/ledger.jsonl` (event `plan_frozen`) before the first stage starts. Changes after that
point are made only as dated amendments appended at the end, never by editing this text.

## 0. What this campaign is, and what the author already knew

A complete rerun of the repaired pipeline in `campaigns/c2026-10-03/`, separate from the historical results in `data/`, `results/` and `figures/`.
The author had seen the HISTORICAL results before writing this plan (surrogate accuracy by split, GA trajectories, the docking result, the H5
experiment). Wherever an operational definition below was shaped by something seen historically, it is marked **[amendment motivated by historical
observation]**, so that the reader can tell a pre-specified endpoint from a convenient one. No campaign output has been inspected.

Frozen inputs (copied into `inputs/`, hashes in `manifests/inputs.json`, each equal to the hash recorded for the historical file):
* `inputs/raw/chembl_CHEMBL301_activities.csv`, sha256 `e42b76cf2427ffa7731dd459b66818d38f8767c591b24280bd21eeb2f18a8011` (the cached ChEMBL download;
  no newer release is downloaded)
* `inputs/structures/4KD1.pdb`, sha256 `f902a65172345488551d1eb14659d2e24094aa82b5b882d88ddecccdf0811f5c`
* `inputs/structures/4kd1_1QK_crystal.sdf`, sha256 `154ef0db46fed220680b8e0943dc270d8416353af89adc3754e0f4373a2e8c08`
Source code: content hash of `src/`, `scripts/`, `tests/` at initialisation `c7f4d920af6fd7c7ef29926afb21693cb176ce49b57b58bfca59244e20b93034`
(78 files; the Git working tree is not clean, so the content hash, not the commit, identifies the code; per-file hashes in
`manifests/source_identity_at_init.json`). Environment captured anew at initialisation (`manifests/environment_at_init.json`): Python 3.11.16, RDKit
2026.3.1, Vina 1.2.7, Meeko 0.8.0, PoseBusters 0.6.5, scikit-learn 1.9.1, NumPy 2.4.6, pandas 3.0.6, macOS arm64, 8 cores. This is a NEW capture; it says
nothing about the environment that produced the historical results. Historical artifacts: `historical_record/historical_hashes_2026-10-03_pre-campaign.json`.

## 1. Original hypotheses (historical wording, preserved verbatim from CLAUDE.md)

> When a molecular optimizer is driven by a machine-learned activity surrogate under multi-objective drug-likeness constraints, does improvement in the
> objective reflect progress toward the target, or exploitation of the surrogate's extrapolation error?

- **H1** Predicted pActivity rises monotonically across GA generations. *(Expected, trivial — the setup, not the finding.)*
- **H2** That rise coincides with increasing distance from the surrogate's training distribution: max Tanimoto to training set *falls* as predicted activity climbs.
- **H3** The gain does not transfer to an orthogonal signal: Vina scores of optimized molecules do not exceed those of held-out known actives.
- **H4** A GA driven by a **label-scrambled** surrogate shows a similar-looking optimization curve — so the curve alone is not evidence of anything.
- **H5** Constraining the GA to stay near the training distribution attenuates H1 but improves transfer in H3. Measure the exchange rate.

## 2. Design

**Preparation policies.** `legacy` scores edited molecules as they are; `corrected` standardizes each candidate like the training set, accepts only closed-shell
structures, deduplicates after standardization, and scores that one representation (cdk2moo/candidates.py). The policy comparison changes several things at
once (standardization, closed-shell rule, post-standardization deduplication). **It is evaluated only as that combined intervention; no component ablation is
run.** Each policy has its own scope; headline results are never pooled across policies.

**Matched initialization.** One eligible-starts list (cdk2moo/starts.py; stage `shared/eligible_starts`) is defined under the stricter policy: inside the size
window, unchanged by standardization, closed-shell, unique after standardization. Both policies, and every constrained run and its unconstrained twin, draw
their starting molecules from it with a seeded generator keyed by (GA seed, allowed set). Paired comparisons are paired by GA seed and start population.

**Surrogates.** Random forest, 300 trees, max_features 0.33, `random_state` 42; scrambled surrogate: one fixed label permutation (seed 20260930); ECFP4
(Morgan radius 2, 2048 bits, binary, no chirality). Trained once per experiment and reused for every GA seed and arm. Experiments: **primary** (all curated
molecules); **molecule-level holdout** (without the 204 top-decile actives); **scaffold-level holdout** (without the 9 scaffold families with at least 5
top-decile members: 260 molecules, 93 of them actives); **five scaffold-split surrogates** (split seeds 42-46; training sets of the Stage 3 scaffold split).

**GA.** Population 100, 100 children per generation, 50 generations, mutation probability 0.5, size window 20-39 heavy atoms, four arms (multi_real,
multi_scrambled, activity_only, druglike_only), GA seeds 42-46, fitness terms scaled 0-1 and combined by geometric mean, activity scaled linearly over
pActivity 4-9. Unchanged from the historical experiments.

**Constrained GA.** Hard floor tau on the raw ECFP4 max-Tanimoto to the experiment's training set, tau in {0.4, 0.5, 0.6, 0.7}, arms multi_real and activity_only,
experiments primary / molecule / scaffold, GA seeds 42-46, each paired with its unconstrained twin.

**Docking.** 4KD1, hydrogens at pH 7.4, 22 A cube on the crystal ligand, Vina exhaustiveness 8, 9 poses, Vina seeds 42/43/44 (43 and 44 only for the replicate
subsets), ligand-preparation seed 42, redocking control first (pass: top pose RMSD <= 2 A). Docked sets are organised by group so that every comparator states what
the driving surrogate saw (cdk2moo/docksets.py):

| group | surrogate that drove the GA | comparator actives | role of the comparator |
|---|---|---|---|
| P | all curated molecules | 100 random top-decile molecules | **in-sample**: the surrogate was trained on them; NOT held out |
| M | without the 204 top-decile actives | 100 random of the 204 held-out actives | **held-out** |
| S | without 9 scaffold families | all 93 held-out actives | **held-out** |

Each group also docks 100 random reference molecules, 100 property-matched decoys, and 100 molecules from the final generation of each of the four arms (group P
also the pooled Pareto front). Reference sets are identical in both scopes and are docked once (shared store). Six molecules from each of five sets per group
are re-docked with two more seeds. **Only groups M and S can assess H3.**

## 3. Operational definitions and decision rules

Unit of replication. For GA-level statistics the unit is the GA RUN (seed): n = 5 per arm and condition, paired by seed. Runs that share a surrogate or overlap in
training data are not independent; the five scaffold-split surrogates share most of their training molecules, so robustness summaries use the SPLIT as the unit
(n = 5, each split summarised by the median over its 5 GA seeds). With n = 5 an exact paired sign-flip test cannot go below two-sided p = 0.0625, so intervals,
ranges and sign counts are the primary uncertainty statements and p-values are secondary. Molecule-level docking statistics use a cluster bootstrap
(2,000 resamples, generator seed 42, 95% percentile intervals) with clusters = GA seed for generated molecules and decoys, and Murcko scaffold for reference
molecules, resampled independently in each of the two samples compared. Surrogate-accuracy results over 10 split seeds are descriptive (overlapping test sets).
**Weighting.** Generated-molecule metrics are reported as DISTINCT-molecule metrics (distinct standardized isomeric structure, represented by its standardized
SMILES) as primary and population-weighted (raw SMILES, every occurrence) as secondary. Only generated molecules (birth_generation > 0) enter statistics about generated
chemistry; generation 0 is the starting reference.

**H1 (assessed: primary experiment, activity-driven arms multi_real and activity_only, each policy).** Per run, m_g = median real-surrogate prediction of generated
molecules in the population at generation g (g = 0: the starting molecules; generations with fewer than 10 generated molecules are excluded from step counts and
reported). Quantities: net rise m_50 - m_0; fraction of consecutive generation pairs with m_(g+1) >= m_g; Spearman(g, m_g). *Rises*: net rise > 0 in all 5 runs and
Spearman >= 0.9 in at least 4. *Monotonically* (the literal word): every consecutive step non-decreasing in every run (fraction = 1); reported separately and stated
as stronger than "rises". Non-activity arms are controls.

**H2 (assessed: primary experiment, activity-driven arms).** Two operationalizations, both reported.
(i) *H2-literal*: per run, Spearman between m_g and s_g over generations 1-50, s_g = median raw max-Tanimoto to the training set of the generated molecules; holds if negative
in at least 4 of 5 runs. (ii) *H2-attributable* **[amendment motivated by historical observation: similarity fell in a drug-likeness-only arm]**: paired final-generation
difference in median raw max-Tanimoto, arm minus druglike_only of the same seed and start; supported if the arm ends farther than the baseline (difference < 0) in at least
4 of 5 seeds. Both are repeated with the size-conditioned percentile (cdk2moo/similarity.py) and with comparison inside heavy-atom strata, because raw similarity depends
on size. H2 as worded is the literal version; the attributable version addresses a different (confound-controlled) question and is labelled as an amendment.

**H3 (assessed ONLY in groups M and S, each policy).** *Validity gate.* An orthogonal signal can only show "no transfer" or "transfer" if it carries activity information. Gate, per
group and policy: the lower bound of the 95% cluster-bootstrap interval of AUC(held-out actives vs matched decoys), score = -Vina, exceeds 0.5 both raw and when stratified by
heavy-atom stratum (20-24, 25-29, 30-34, 35-39; stratum-pair-weighted mean AUC). A pass is labelled "weak" if the point AUC < 0.65 and "moderate" otherwise; a pass does not
establish usefulness. *Endpoint if the gate passes:* P_sup(arm) = probability that a random distinct generated molecule of the arm scores better (more negative) than a random
held-out active, cluster-bootstrap interval, raw and stratified. For an arm: H3 "not exceeded" if the interval's lower bound <= 0.5; "exceeded" if it is > 0.5. *If the gate
fails:* H3 is reported as **not assessable with Vina on this chemistry**; P_sup is still tabulated, labelled descriptive and not interpretable as transfer. Group P is reported as
an in-sample comparison (exploratory), never as H3. Secondary docking analyses: score vs measured pActivity (reference sets), size dependence (slope over all molecules and
within strata), replicate noise (3 seeds), PoseBusters (three states), hinge proximity (distance-only: ligand N or O within 3.5 A of Glu81 O or Leu83 N/O).

**H4 (assessed: primary experiment, multi_real vs multi_scrambled, each policy).** Per run, F_g = median prediction of the surrogate the arm FOLLOWS among generated
molecules (real surrogate for multi_real; scrambled surrogate for multi_scrambled). "Similar-looking curve" is operationalized as Spearman(g, F_g) >= 0.9 in at least 4 of 5
runs for both arms. Magnitude is reported separately: net rise of F, and the rise standardized by the standard deviation of that surrogate's predictions over its own training
molecules, with the ratio scrambled / real (no threshold). The interpretive clause ("so the curve alone is not evidence") is addressed by the non-transfer quantity: the real-surrogate
prediction of multi_scrambled molecules minus that of druglike_only (paired by seed), compared with the same quantity for multi_real. H4 is a statement about curves, not about
surrogate quality.

**H5 (assessed: all three experiments, arms multi_real and activity_only, each policy).** *Constraint*: the hard floor above (one operationalization of "stay near the training
distribution"). *Attenuation of H1*: paired difference in final median predicted activity, constrained minus unconstrained twin; attenuation if negative in at least 4 of 5 seeds,
shown across the four floors. *Transfer in H3 (original wording)*: requires Vina to be validated; **pre-registered contingency**: constrained molecules (tau = 0.6, arms multi_real and
activity_only, groups M and S, 100 distinct final-generation molecules each, same selection rules as other sets) are docked ONLY IF the H3 validity gate passes for that group and
policy; otherwise H5 transfer in H3 terms is reported as not assessable. *Amended exploratory transfer* **[amendment motivated by historical work]**: share of distinct final-generation
molecules at ECFP4 proximity >= 0.6 to a held-out active (a fingerprint-closeness measure, not activity), reported in all cases and labelled as such. *Exchange rate*: defined only
where attenuation exists: change in the transfer endpoint per unit of predicted activity given up, per floor, seed-level bootstrap interval; otherwise "undefined" with the reason.

## 4. Analysis tiers

**Primary.** H1, H2 (both versions), H4, H5 attenuation, the H3 validity gate, and (if the gate passes) H3 in groups M and S: per policy, primary experiment or holdout experiments as
stated. **Secondary.** Five scaffold-split robustness (split as unit; a claim is "stable" if it holds in the primary run and at least 4 of 5 splits); applicability-domain reliability
(Stage 4 out-of-sample error vs raw similarity and vs size percentile); generated molecules placed on that curve; structural alerts (PAINS, Brenk; recorded, not filtered);
validity diagnostics (radical counts, rejection reasons, tautomer-enumeration status); policy comparison (side by side); Pareto analysis; surrogate evaluation across 10 split
seeds with scrambled, mean and size-only controls. **Exploratory.** Anything not listed (hinge proximity patterns, chemotype/diversity descriptions, scaffold-family-level
observations, in-sample group P docking). Exploratory observations are labelled as such in the report and are not used to support hypothesis-level statements.

## 5. Missing data and failures
A docking job that fails to prepare or dock stays in the ledger and in every denominator as a failure (counts reported per set); score statistics use completed jobs and state n. A pose-export
failure is a failure. PoseBusters results have three states (pass, fail, unevaluable) with the denominator stated; unevaluable poses are never counted as passes, and primary Vina analyses are
repeated excluding poses that fail or are unevaluable. A molecule whose standardization fails is counted separately and kept as its own identity. A GA run that errors is recorded as failed;
the stage is fixed (if an implementation fault) and rerun in full with the same seeds, and the failed attempt is kept under `superseded/`.

## 6. Experiments, budgets and stopping conditions
Stages (67, see `src/cdk2moo/campaign_driver.py`): shared curation, fingerprints and splits, eligible starts, surrogate evaluation (10 split seeds x 3 split types x 4 models), holdout setup,
receptor preparation and redocking; then per policy: primary GA (4 arms x 5 seeds), molecule- and scaffold-level holdout GA, five scaffold-split GA, robustness and positive-control
analyses, Pareto, applicability domain, alerts, constrained GA (3 experiments x 2 arms x 4 floors x 5 seeds), docking-set selection, docking (about 2,000 distinct molecules per policy
plus replicates, shared store), docking analysis, PoseBusters, hinge, figures, Table 1, audits. Budget: unlimited by the user's authorisation. Seeds, parameters and thresholds are frozen as
above. Stopping: every stage complete or recorded as blocked with its exact blocker; there is no early stopping on the appearance of plots, and no parameter changes in response to results.

## 7. What each experiment can and cannot assess

| experiment | can assess | cannot assess |
|---|---|---|
| primary GA, 4 arms, both policies | H1, H2, H4; effect of the combined policy change | whether the optimized molecules bind CDK2 (no experimental measurement exists) |
| holdout GAs (M, S) | transfer to held-out actives by fingerprint proximity; the H3 comparison when docked | whether proximity implies activity; independence of siblings left in training |
| constrained GAs | H5 attenuation; amended proximity transfer; exchange rate where defined | H5 transfer in H3 terms unless the docking contingency is triggered |
| docking (groups M, S) | whether Vina discriminates held-out actives (gate) and, if so, H3 | binding affinity; Vina is a rigid-receptor score without rigorous entropy |
| docking (group P) | whether Vina tracks activity for in-sample actives | H3 (the comparator is not held out) |
| five scaffold-split surrogates | stability of primary claims to the surrogate's training subset | independence (training sets overlap) |
| leakage audits | train/test and control overlap | temporal generalisation (no temporal split exists; `primary_document` is not a verified earliest publication) |

## 8. Known limitations accepted in advance
Assay heterogeneity (ATP concentration mostly unrecorded; many assays on CDK2/cyclin A while docking uses monomeric CDK2); no measurable noise floor; the "held-out" controls leave close
relatives in training (88% and 76% of held-out actives had a remaining-training neighbour at Tanimoto >= 0.6 historically; recomputed in this campaign); ECFP4 cannot tell molecules with
identical fingerprints apart and similarity-based models smooth over activity cliffs; GA edits are chemically naive; QED rewards resemblance to existing drugs; tautomer enumeration is
incomplete for about 1% of molecules; RDKit's tautomer step strips some stereochemistry; nothing is validated experimentally.

## 9. Amendments
(none yet)

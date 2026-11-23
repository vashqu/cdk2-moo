# Campaign c2026-11-22: execution and analysis plan (frozen before any campaign stage runs)

This file is hashed (`plan/plan.sha256.json`) before the first stage starts; `campaign_run.py run` refuses to start without a matching hash. Changes after that
point are made only as dated amendments in `CAMPAIGN_PLAN_AMENDMENTS.md`, never by editing this text.

## 0. What this campaign is, and what the author already knew

A second campaign, run after an independent review of campaign `c2026-10-03` found analysis and code problems (`docs/known_issues.md`). It repeats the experiment with those
problems repaired and with a stronger design. **The author has seen every result of `c2026-10-03`** (trajectories, the failed docking validity gate, the H5 table) and the
review's findings; the plan below was written with that knowledge. Wherever a definition below was shaped by something seen earlier it is marked **[motivated by c2026-10-03]**.
No output of this campaign has been inspected.

**Changes from c2026-10-03, and why**

| Change | Reason |
|---|---|
| Only the `corrected` preparation policy is run | `legacy` scored candidates as edited (radical-bearing and duplicate candidates). Its results stay in `c2026-10-03` as the historical comparison. The bundled legacy-vs-corrected question is **not** re-asked. |
| 10 GA seeds per arm (42-51) instead of 5 | With n = 5 the smallest exact paired p-value was 0.0625 and every verdict rested on "4 of 5 runs". |
| One label permutation **per GA seed** for the scrambled control (seed 42 keeps the permutation of c2026-10-03; permutation seed = 20260930 + GA seed - 42) | One fixed permutation made the control vary only the search. Now it also samples over permutations. |
| Decision rules read "at least 80% of runs" (8 of 10) where they read "4 of 5" | Same proportion, scaled to the new n. |
| Primary counting unit: distinct standardized molecules, enforced in code (`cdk2moo/distinct.py`) | c2026-10-03 summarized raw population rows. Under `corrected`, rows and distinct standardized molecules coincide by construction; the check makes that a refusal, not an assumption. |
| Docked poses must keep every stereo element the input specifies (element-wise check at export) | The earlier rule compared full InChIKeys only when every centre was specified; it missed a flipped double bond in 4 of 3,650 saved poses (macrocycles). A failing pose is recorded as a `pose_export_failed` job and stays in the denominator. |
| H5 exchange rate: pairs with missing values are dropped (never zero-filled); the interval is withheld unless the ratio is defined in at least 95% of resamples | The earlier estimator filled missing proximity with zero and conditioned its interval on the denominator's sign. |
| Constrained-GA tables keep all provenance columns | Earlier tables dropped proposal SMILES, parents, origin and tautomer status. |
| Resume checks consumed-input hashes and a hash of `src/` + `scripts/`; the run needs committed code and this frozen plan | Earlier resume checked declared output hashes only; the recorded commit did not contain the code. |

Frozen inputs (copied into `inputs/`, hashes in `manifests/inputs.json`, each equal to the hash recorded for the historical file):
* `inputs/raw/chembl_CHEMBL301_activities.csv`, sha256 `e42b76cf2427ffa7731dd459b66818d38f8767c591b24280bd21eeb2f18a8011` (the cached ChEMBL download; no newer release)
* `inputs/structures/4KD1.pdb`, sha256 `f902a65172345488551d1eb14659d2e24094aa82b5b882d88ddecccdf0811f5c`
* `inputs/structures/4kd1_1QK_crystal.sdf`, sha256 `154ef0db46fed220680b8e0943dc270d8416353af89adc3754e0f4373a2e8c08`

Source code: identified by the git commit and the content hash of `src/` and `scripts/` recorded at initialisation (`manifests/source_identity_at_init.json`); the run is refused
if `src/`, `scripts/`, `tests/`, `pyproject.toml` or `environment.yml` are not committed. Environment: captured at initialisation (`manifests/environment_at_init.json`).

## 1. Original hypotheses (historical wording, preserved)

> When a molecular optimizer is driven by a machine-learned activity surrogate under multi-objective drug-likeness constraints, does improvement in the
> objective reflect progress toward the target, or exploitation of the surrogate's extrapolation error?

- **H1** Predicted pActivity rises monotonically across GA generations. *(Expected, trivial: the setup, not the finding.)*
- **H2** That rise coincides with increasing distance from the surrogate's training distribution: max Tanimoto to training set *falls* as predicted activity climbs.
- **H3** The gain does not transfer to an orthogonal signal: Vina scores of optimized molecules do not exceed those of held-out known actives.
- **H4** A GA driven by a **label-scrambled** surrogate shows a similar-looking optimization curve, so the curve alone is not evidence of anything.
- **H5** Constraining the GA to stay near the training distribution attenuates H1 but improves transfer in H3. Measure the exchange rate.

## 2. Design

**Policy.** `corrected` only: each candidate is standardized like the training set (cleanup, largest fragment, neutralise, canonical tautomer, organic-element check), accepted only
if closed-shell, deduplicated after standardization, and scored in that one representation (`cdk2moo/candidates.py`).

**Matched initialization.** One eligible-starts list (`shared/eligible_starts`: inside the size window, unchanged by standardization, closed-shell, unique) supplies the starting
molecules of every run, with a seeded generator keyed by (GA seed, allowed set). Constrained runs and their unconstrained twins share starts and are paired by GA seed.

**Surrogates.** Random forest, 300 trees, max_features 0.33, `random_state` 42; ECFP4 (Morgan radius 2, 2048 bits, binary, no chirality). The real surrogate is trained once per
experiment and reused for every seed and arm. The scrambled surrogate of GA seed *s* is trained on a label permutation with seed 20260930 + (s - 42) and a forest with the same seed.
Experiments: **primary** (all curated molecules); **molecule-level holdout** (without the 204 top-decile actives); **scaffold-level holdout** (without the 9 scaffold families with at least
5 top-decile members); **five scaffold-split surrogates** (split seeds 42-46).

**GA.** Population 100, 100 children per generation, 50 generations, mutation probability 0.5, size window 20-39 heavy atoms, four arms (multi_real, multi_scrambled, activity_only,
druglike_only), **GA seeds 42-51**, fitness terms scaled 0-1 and combined by geometric mean, activity scaled linearly over pActivity 4-9. Unchanged from c2026-10-03 except the seeds.

**Constrained GA.** Hard floor tau on the raw ECFP4 max-Tanimoto to the experiment's training set, tau in {0.4, 0.5, 0.6, 0.7}, arms multi_real and activity_only, experiments
primary / molecule / scaffold, GA seeds 42-51, each paired with its unconstrained twin.

**Docking.** 4KD1, hydrogens at pH 7.4, 22 A cube on the crystal ligand, Vina exhaustiveness 8, 9 poses, Vina seeds 42/43/44 (43 and 44 only for replicate subsets), ligand-preparation
seed 42, redocking control first (pass: top pose RMSD <= 2 A). Groups as in c2026-10-03:

| group | surrogate that drove the GA | comparator actives | role of the comparator |
|---|---|---|---|
| P | all curated molecules | 100 random top-decile molecules | **in-sample** (the surrogate was trained on them; NOT held out) |
| M | without the 204 top-decile actives | 100 random of the 204 held-out actives | **held-out** |
| S | without 9 scaffold families | all 93 held-out actives | **held-out** |

Each group also docks 100 random reference molecules, 100 property-matched decoys, and 100 molecules from the final generation of each of the four arms (group P also the pooled Pareto
front), drawn from all 10 seeds. Six molecules from each of five sets per group are re-docked with two more seeds. **Only groups M and S can assess H3.**

## 3. Operational definitions and decision rules

Unit of replication: the GA RUN (seed); n = 10 per arm and condition, paired by seed. "At least 80% of runs" means at least 8 of 10 (`config.runs_needed`). With n = 10 an exact paired
sign-flip test reaches two-sided p = 0.002; intervals, ranges and run counts remain the primary uncertainty statements, p-values secondary. Runs that share a surrogate or overlap in training data
are not independent: robustness summaries over the five scaffold-split surrogates use the SPLIT as the unit (n = 5, each split summarised by the median over its 10 GA seeds). Docking statistics use
a cluster bootstrap (2,000 resamples, generator seed 42, 95% percentile intervals), clusters = GA seed for generated molecules and decoys, Murcko scaffold for reference molecules. Surrogate-accuracy
results over 10 split seeds are descriptive.

**Counting unit.** Distinct standardized molecules, represented by their stored (standardized, canonical) SMILES; the code refuses populations that are not entirely `corrected` or that repeat a SMILES
within one population. Only generated molecules (birth_generation > 0) enter statistics about generated chemistry; generation 0 is the starting reference (training-set molecules).

**H1 (primary experiment, arms multi_real and activity_only).** Per run, m_g = median real-surrogate prediction of generated molecules in the population at generation g (g = 0: starting molecules;
generations with fewer than 10 generated molecules are excluded from step counts and reported). Net rise m_50 - m_0; fraction of consecutive pairs with m_(g+1) >= m_g; Spearman(g, m_g).
*Rises*: net rise > 0 in all runs and Spearman >= 0.9 in at least 80% of runs. *Monotonically* (the literal word): every consecutive step non-decreasing in every run; reported separately and stated as stronger.

**H2 (primary experiment, activity-driven arms).** (i) *Literal*: per run, Spearman between m_g and s_g over generations 1-50 (s_g = median raw max-Tanimoto to the training set); holds if negative in at
least 80% of runs. (ii) *Attributable* **[motivated by c2026-10-03: similarity fell in a drug-likeness-only arm]**: paired final-generation difference in median raw max-Tanimoto, arm minus druglike_only of the
same seed; supported if negative in at least 80% of runs. Both repeated with the size-conditioned percentile and within heavy-atom strata, because raw similarity depends on size. **[Motivated by
c2026-10-03]** the within-stratum repetition was planned there but not implemented; here it is a required output of stage 13, and its absence is reported as a failure of the stage.

**H3 (groups M and S).** *Validity gate*, per group: the lower bound of the 95% cluster-bootstrap interval of AUC(held-out actives vs matched decoys), score = -Vina, exceeds 0.5 both raw and stratified by
heavy-atom stratum (20-24, 25-29, 30-34, 35-39). A pass is "weak" if point AUC < 0.65, otherwise "moderate"; a pass does not establish usefulness. *If the gate passes*: P_sup(arm) with cluster-bootstrap
interval, raw and stratified; "exceeded" if the lower bound > 0.5, otherwise **"superiority not established"** (this is not "no difference"). *If the gate fails*: H3 is **not assessable with Vina on this
chemistry**; P_sup is tabulated as descriptive only. In c2026-10-03 the gate failed in every group; the same outcome is plausible here and is an acceptable result. Group P is an in-sample comparison, never H3.
Secondary: score vs measured pActivity, size dependence, replicate noise, PoseBusters (three states), hinge proximity (distance only), and the count of poses rejected by the stereo check.

**H4 (primary experiment, multi_real vs multi_scrambled).** Per run, F_g = median prediction of the surrogate the arm FOLLOWS (real for multi_real; that seed's scrambled surrogate for multi_scrambled).
"Similar-looking curve": Spearman(g, F_g) >= 0.9 in at least 80% of runs for both arms. Magnitude separately: net rise of F, and the rise divided by the standard deviation of that surrogate's predictions over
its own training molecules (per seed's own scrambled surrogate), with the ratio scrambled / real (no threshold). Non-transfer quantity: real-surrogate prediction of multi_scrambled molecules minus
druglike_only (paired by seed), compared with multi_real's.

**H5 (all three experiments, arms multi_real and activity_only).** *Attenuation*: paired difference in final median predicted activity, constrained minus unconstrained twin; attenuation if negative in at least 80% of
seeds, shown per floor. *Transfer in H3 terms*: constrained molecules (tau = 0.6, both arms, groups M and S, 100 distinct final-generation molecules each) are docked ONLY IF the H3 validity gate passes for that group;
otherwise not assessable. *Amended exploratory transfer* **[motivated by c2026-10-03]**: share of final-generation molecules at ECFP4 proximity >= 0.6 to a held-out active (fingerprint closeness, not activity), reported
and labelled as such; it shares geometry with the imposed floor and is not an independent endpoint. *Exchange rate*: computed only where attenuation holds by the rule above; per-seed pairs with a missing value are
dropped (never zero-filled); the bootstrap interval is withheld unless the ratio is defined in at least 95% of resamples; otherwise reported "undefined" with the reason. **No headline claim rests on the exchange rate.**

## 4. Analysis tiers

**Primary.** H1, H2 (both versions), H4, H5 attenuation, the H3 validity gate and, if the gate passes, H3 in groups M and S. **Secondary.** Five scaffold-split robustness (split as unit; "stable" if it holds in the primary
run and at least 4 of 5 splits); applicability-domain reliability; structural alerts (recorded, not filtered); validity diagnostics (radical counts, rejection reasons, tautomer-enumeration status, and the share of stored
SMILES that change when standardized again); Pareto analysis; surrogate evaluation over 10 split seeds with scrambled, mean and size-only controls. **Exploratory.** Anything not listed (hinge patterns, chemotype
descriptions, group P docking). Exploratory observations are labelled as such and not used for hypothesis-level statements.

## 5. Missing data and failures
A docking job that fails to prepare, dock or export a stereo-faithful pose stays in the ledger and in every denominator as a failure (counts reported per set); score statistics use completed jobs and state n.
PoseBusters results have three states (pass, fail, unevaluable); unevaluable poses are never passes; primary Vina analyses are repeated excluding failing or unevaluable poses. A molecule whose standardization fails is
counted separately. A GA run that errors is recorded as failed; an implementation fault is fixed and the stage rerun in full with the same seeds, the failed attempt kept under `superseded/`.

## 6. Budget and stopping
Stages: see `src/cdk2moo/campaign_driver.py` (`corrected` policy: shared stages, then GA, holdout, robustness, constrained, docking-set selection, docking, analysis, PoseBusters, hinge, figures, audits). Expected wall time
about a day on 8 cores (GA stages scale with the 10 seeds; docking about 3-4 h). Seeds, parameters and thresholds are frozen above. Stopping: every stage complete or recorded as blocked with its exact blocker; no early
stopping on the appearance of plots; no parameter changes in response to results. A stage failure is fixed only if it is an implementation fault, with the fix committed before the stage is rerun.

## 7. What each experiment can and cannot assess

| experiment | can assess | cannot assess |
|---|---|---|
| primary GA, 4 arms | H1, H2, H4 | whether the optimized molecules bind CDK2 (no experimental measurement exists) |
| holdout GAs (M, S) | transfer to held-out actives by fingerprint proximity; the H3 comparison when docked | whether proximity implies activity; independence of siblings left in training |
| constrained GAs | H5 attenuation; amended proximity transfer; exchange rate where defined | H5 transfer in H3 terms unless the docking contingency is triggered |
| docking (groups M, S) | whether Vina discriminates held-out actives (gate) and, if so, H3 | binding affinity |
| docking (group P) | whether Vina tracks activity for in-sample actives | H3 (comparator not held out) |
| five scaffold-split surrogates | stability of primary claims to the training subset | independence (training sets overlap) |
| leakage audits | train/test and control overlap | temporal generalisation |

## 8. Known limitations accepted in advance
Assay heterogeneity (ATP concentration mostly unrecorded; many assays on CDK2/cyclin A while docking uses monomeric CDK2); censored `>` records dropped; no measurable noise floor; "held-out" controls leave close
relatives in training; ECFP4 cannot tell molecules with identical fingerprints apart and similarity-based models smooth over activity cliffs; GA edits are chemically naive; QED rewards resemblance to existing drugs;
tautomer enumeration is incomplete for about 1% of molecules (so re-standardizing a stored SMILES can occasionally change it); the corrected-versus-legacy comparison is not repeated here; nothing is validated experimentally.

## 9. Amendments
(none yet)

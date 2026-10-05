# CDK2 multi-objective molecular optimization

**Question.** When a molecular optimizer is driven by a machine-learned activity surrogate under drug-likeness constraints, does improvement in the objective reflect progress
toward the target, or exploitation of the surrogate's extrapolation error?

A negative or null result is a valid outcome here. Judgment calls are in [`decisions.md`](decisions.md) (append-only). Project rules and known confounds are in
[`CLAUDE.md`](CLAUDE.md). The current evidence is campaign **`c2026-10-03`** (`campaigns/c2026-10-03/`); the audit hand-off with the claim-to-evidence table is
[`campaigns/c2026-10-03/HANDOFF.md`](campaigns/c2026-10-03/HANDOFF.md).

## Where things are (historical vs current: never mixed)

| What | Where | Status |
|---|---|---|
| Historical results (before 2026-10-03): data, results, figures, docking poses | `data/`, `results/`, `figures/` | **Unchanged.** 132 files hashed before the campaign (`historical_record/`) and re-verified after: none missing; the only changes are the intended rewrites of README.md and CLAUDE.md and the append-only decisions.md. Readable only with `CDK2_CAMPAIGN=historical`. |
| README and CLAUDE.md as they stood before the campaign | `historical_record/*_before_campaign_c2026-10-03.md` | Preserved; hashes equal the recorded ones. `verify-historical` now reports exactly these two and `decisions.md` (append-only) as changed. |
| Current campaign: frozen inputs, plan, ledger, manifests | `campaigns/c2026-10-03/{inputs,plan,manifests}/` | Plan hashed before any stage ran; amendments in `plan/CAMPAIGN_PLAN_AMENDMENTS.md`. |
| Shared stages (curation, splits, surrogate evaluation, holdout setup, receptor, redocking) | `campaigns/c2026-10-03/shared/` | Policy-independent. |
| Results under the **legacy** preparation policy (candidates scored as edited) | `campaigns/c2026-10-03/legacy/{results,figures}/` | Headline results are never pooled across policies. |
| Results under the **corrected** policy (standardize like training data, closed-shell only, deduplicate after standardization) | `campaigns/c2026-10-03/corrected/{results,figures}/` | The two policies are compared only as that combined intervention. |
| Docking job store (shared; a molecule docked once is reused) | `campaigns/c2026-10-03/docking_store/` | 3,660 job records (3,654 ok, 6 recorded `prep_failed`), poses hashed. |
| Failed or superseded attempts | `campaigns/c2026-10-03/superseded/`, ledger | Kept as evidence, never reused. |

## Reproduce

From the repo root, in the existing `cdk2` conda environment (no package changes; Python 3.11, RDKit 2026.3.1, Vina 1.2.7, Meeko 0.8.0, PoseBusters 0.6.5; the environment capture for this
campaign is `manifests/environment_at_init.json` and says nothing about the environment of the historical results):

```
python -m unittest discover -s tests                                   # 159 tests, standard library only
python scripts/campaign_run.py init   --campaign c2026-10-03          # freeze inputs (hash-checked against the historical files), write the ledger
python scripts/campaign_run.py run    --campaign c2026-10-03 --cores 8 # all 69 stages in dependency order; resumable; docking was about 8.5 h on 8 cores (5.2 h legacy, 3.3 h corrected after reuse)
python scripts/campaign_run.py status --campaign c2026-10-03          # ledger view
python scripts/campaign_run.py verify-historical                       # historical files unchanged?
```

A single stage runs as `CDK2_CAMPAIGN=c2026-10-03 CDK2_SCOPE=legacy python scripts/<stage>.py`. A script run without a campaign stops; it never falls back to historical files.
The campaign was executed in several `run --only <patterns>` invocations (GA stages, analysis stages, docking stages, figures); the ledger has every attempt, including eleven failed attempts
(holdout_eval x4, dock x1, table1 x2, figures x4) that were repaired (D-42) and re-run. `CAMPAIGN_PLAN_AMENDMENTS.md` lists which code changed after the plan was frozen (A7).

## What was run

Both policies, same frozen inputs (the cached ChEMBL download; no newer release): 2,016 curated molecules; random-forest surrogate evaluated on random, scaffold and paper splits (10 seeds each) with scrambled-label,
mean and size-only controls; GA with four arms (multi_real, multi_scrambled, activity_only, druglike_only), 5 GA seeds, for three surrogates (trained on all data / without 204 held-out actives /
without 9 scaffold families) plus five scaffold-split surrogates; constrained GA (similarity floor 0.4-0.7, arms multi_real and activity_only, three experiments); Pareto, applicability-domain, structural-alert and
validity diagnostics; docking of 2,033 (legacy) and 1,999 (corrected) distinct molecules into 4KD1 with Vina (+2 more seeds for 89-94 replicate molecules), PoseBusters (22 checks), hinge proximity,
leakage audits. Matched initialization: both policies and each constrained run draw their starting molecules from one eligible-starts list.

Not run, and why: the **docking-based H5 contingency** (constrained molecules docked) did not trigger because the pre-registered H3 validity gate failed (below); a **random-ChEMBL** reference sample (the
"random" sets are random training molecules, a stand-in); component ablation of the three changes bundled in the corrected policy; slide deck (author's task).

## Results (campaign `c2026-10-03`; every number is per policy, never pooled)

Unit of replication is the GA run (n = 5 per arm, paired by seed). With n = 5 the smallest possible exact two-sided sign-flip p is 0.0625, so ranges and counts of runs are the primary evidence.
"Legacy / corrected" below always means the two policy scopes.

**Surrogate (shared).** Test R2 median over 10 split seeds, random / scaffold / paper: RF 0.64 / 0.49 / 0.30 (Spearman 0.81 / 0.71 / 0.54); label-scrambled -0.12 / -0.13 / -0.08; size-only 0.05 / 0.01 / -0.03; mean -0.01.
The surrogate is an instrument; a scaffold-split R2 of about 0.5 says it is useful near known chemistry and weaker away from it.

| Hypothesis (as worded in the project) | Verdict | Evidence (legacy / corrected) |
|---|---|---|
| **H1** predicted activity rises monotonically | **Rises: yes. Monotonic: only activity_only.** | Net rise of median real prediction, multi_real +1.07 / +1.13, activity_only +1.91 / +1.93 (5/5 runs, Spearman >= 0.9 in 5/5). Every step non-decreasing: activity_only 5/5 runs; multi_real 0/5 (86% of steps non-decreasing). |
| **H2** the rise coincides with falling similarity to training | **Holds literally for multi_real; not attributable to the activity objective; reversed for activity_only** | multi_real: similarity falls with prediction in 4/5 runs (median Spearman -0.95 / -0.92), but druglike_only (no activity term) ends at the same similarity (0.36); multi_real minus druglike_only paired final similarity +0.024 / +0.007, farther in 1/5 and 2/5 seeds. activity_only: similarity rises to ~0.70, +0.35 above baseline, negative in 0/5 runs. |
| **H3** gain does not transfer to Vina-scored held-out actives | **Not assessable.** | Validity gate failed in every group and policy: AUC(held-out actives vs matched decoys, -Vina) raw [95% cluster-bootstrap] M 0.59 [0.46, 0.72] / 0.55 [0.43, 0.69]; S 0.57 [0.43, 0.77] / 0.54 [0.39, 0.77]; in-sample P 0.57 [0.46, 0.68] / 0.60 [0.49, 0.70]. Vina vs measured pActivity Spearman 0.01-0.10, intervals include 0. Docking score tracks heavy-atom count (-0.04 to -0.07 kcal/mol per atom). |
| **H4** scrambled-label GA looks like it optimizes | **Holds for the curve; no transfer** | Followed-prediction Spearman >= 0.9 in 5/5 runs for both arms; scrambled rise is smaller (standardized 0.75 vs 1.06; 0.77 vs 1.12). The real surrogate's prediction for multi_scrambled molecules minus druglike_only: -0.02 / -0.13 (above baseline in 2/5 seeds), against +1.56 / +1.64 for multi_real. |
| **H5** a similarity floor attenuates H1 and improves H3 transfer | **Attenuation: partial, arm-specific. Transfer: not assessable in H3 terms.** | multi_real: the median final prediction is at or above the unconstrained twin at every floor (-0.01 to +0.46), never attenuated in >= 4/5 seeds. activity_only: attenuation in >= 4/5 seeds only at some floors (floor 0.7 in the molecule and scaffold experiments of both policies and the primary experiment of the corrected policy; median differences -0.04 to -0.37 pActivity). Exchange rate (change in proximity share per pActivity given up) has an interval excluding 0 only for scaffold-holdout activity_only at floor 0.7 in both policies (legacy 39.9 [3.5, 86.9], corrected 34.8 [8.3, 94.5] pp per pActivity unit) and is otherwise undefined or spans 0. The amended exploratory transfer measure is fingerprint proximity to held-out actives, not activity; by construction it rises when molecules are forced near training analogues of the held-out actives (88% / 76% of held-out actives have such a neighbour). |

**Scaffold-split robustness** (split as the unit, n = 5): four of six pre-declared claims are stable in both policies; "multi_real keeps raw similarity >= druglike_only" (3/5 splits) and "multi_real and
activity_only keep fewer distinct scaffolds" (2/5 legacy, 3/5 corrected) are fragile. This analysis uses population-weighted medians (the older definition).

**What the activity pressure does to the molecules (descriptive).** multi_real converges on the smallest molecules the window allows (median 21 heavy atoms, window minimum 20) with QED 0.91, 39 distinct scaffolds per 100;
activity_only grows large (median 35 atoms, MW ~490, QED ~0.3) and stays near known chemistry (similarity ~0.70). Because docking scores fall with size, the two arms occupy opposite ends of the size confound.

**Preparation policy.** The two policies reach the same verdict on every row above. They differ in validity: under the legacy policy 149 of 12,095 distinct generated (arm, SMILES) molecules (activity_only 106) carry an unpaired-electron (radical) atom,
and 8 docked legacy poses fail the PoseBusters `no_radicals` check; under the corrected policy there are none, by construction. The combined intervention is not decomposed.

**Held-out proximity (fingerprint closeness, not activity).** Share of distinct generated molecules at ECFP4 >= 0.6 to a held-out active, molecule-level holdout, legacy / corrected: activity_only 84% / 90%, multi_real 2% / 2%, druglike_only 0% / 0%,
starting molecules 28%; remaining training molecules alone: 26%. Held-out actives with an activity-only neighbour at that closeness: 107 of 204 / 104 of 204. At the scaffold-level holdout: activity_only 16% / 7%.
Siblings left in training already give 88% (molecule) and 76% (scaffold) of held-out actives a training neighbour at >= 0.6.

**Docking checks.** Redocking of dinaciclib: top-pose RMSD 0.63-0.65 A over three Vina seeds (control passed, < 2 A), so docking is mechanically sound; that says nothing about activity. Engine noise (89 molecules, 3 seeds): median SD 0.012 kcal/mol,
maximum range 2.63. PoseBusters: 1,991 of 2,032 legacy and 1,967 of 1,996 corrected audited poses pass all 22 checks (a docking engine avoids clashes by construction; weak evidence). Hinge proximity (distance only: ligand N/O within 3.5 A of Glu81 O or Leu83 N/O)
is 76-92% for activity_only (by group and policy) and 88-95% for held-out actives, against 38-55% for druglike_only (exploratory; not a verified hydrogen bond, and not used for any hypothesis).

**Not answered here.** Whether any generated molecule is active against CDK2 (nothing was measured); whether the "held-out" actives are genuinely unseen chemistry (most have training siblings); H3 and the H5 transfer claim (no validated orthogonal signal);
which of the three changes in the corrected policy matters; temporal generalization.

Figures and tables per policy: `campaigns/c2026-10-03/<policy>/figures/final/` (fig1 surrogate, fig2 trajectories, fig3 Pareto, fig4 similarity vs prediction with the surrogate's measured error, fig5 docking per group P/M/S, fig6 pose checks, `table1.md`).
Every title states campaign, policy, sample sizes, weighting (distinct molecules unless stated) and uncertainty.

## Limitations

Assay heterogeneity (ATP concentration mostly unrecorded; many assays on CDK2/cyclin A while docking uses monomeric CDK2); no measurable noise floor; censored ">" records dropped; Vina is a rigid-receptor score
with no rigorous entropy and carries no activity signal on this chemistry; ECFP4 cannot tell molecules with identical fingerprints (e.g. stereoisomers) apart and similarity-based models tend to smooth over
activity cliffs; GA edits are chemically naive; QED rewards resemblance to existing drugs; "held-out" controls leave close relatives in training; tautomer enumeration is incomplete for ~1% of molecules and RDKit's tautomer step
strips some stereochemistry; n = 5 GA runs per arm; docking-set sizes (100 per set) are modest; seeds, parameters and thresholds were frozen before the campaign. Nothing is validated experimentally.
Implementation correctness (159 tests, hash-verified stores) is not scientific validity: the tests establish that metrics are computed as defined and that stores are consistent, not that the experiment answers its question.

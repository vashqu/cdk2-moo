# CDK2 multi-objective molecular optimization

**Question.** When a molecular optimizer is driven by a machine-learned activity surrogate under drug-likeness constraints, does improvement in the objective reflect progress toward the target, or exploitation of the surrogate's extrapolation error?

**Current evidence: campaign `c2026-10-03`** (`campaigns/c2026-10-03/`), documented as of 2026-10-05 after an independent review (`../review-2026-10-05/REVIEW.md`). The project is a methodological case study. It does **not** establish whether any generated molecule is active against CDK2, and it did not measure the
surrogate's error on generated molecules. Several analysis problems found in review are **unresolved**; see [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md). Affected quantities are marked provisional below. This README was rewritten in a documentation-only pass: no code, test, analysis or figure was run or changed.

Start with [`presentation/START_HERE.md`](presentation/START_HERE.md). Claims with their sources: [`presentation/CLAIMS_TO_EVIDENCE.md`](presentation/CLAIMS_TO_EVIDENCE.md). Figures with captions and limitations: [`presentation/FIGURE_INDEX.md`](presentation/FIGURE_INDEX.md).
Judgment calls: [`decisions.md`](decisions.md) (append-only). Project rules and the original hypotheses: [`CLAUDE.md`](CLAUDE.md).

## Methods in brief

* **Data.** CHEMBL301 (human CDK2) IC50 records from a cached ChEMBL download (5,367 rows; no newer release used), curated to 2,016 standardized structures: exact-relation `=` records only (values reported as ">" were dropped, a modelling choice that biases the set toward actives), mutant-protein assays removed, median pActivity over each structure's **records** (records are not demonstrably independent experiments).
* **Surrogate.** ECFP4 fingerprints (radius 2, 2048 bits) with a random forest (300 trees); evaluated on random, scaffold and "primary-document" splits, 10 seeds each, with a label-scrambled model, a mean predictor and a size-only regression as controls. The surrogate is an instrument, not a result.
* **Optimizer.** A graph-based genetic algorithm (population 100, 50 generations, heavy-atom window 20-39) with four arms: `multi_real` (activity x QED x SA), `multi_scrambled` (same, driven by a label-scrambled surrogate), `activity_only`, `druglike_only` (QED x SA; no activity term). Five GA seeds (42-46) per arm, for three surrogates (all data; without 204 held-out top-decile actives; without 9 scaffold families), five scaffold-split surrogates, and constrained runs with a hard similarity floor (0.4-0.7).
* **Two preparation policies.** `legacy` scores candidates as edited; `corrected` standardizes like the training data, accepts only closed-shell structures and deduplicates after standardization. The comparison **bundles those three changes**; components were not isolated. Results are reported per policy and never pooled.
* **Docking.** AutoDock Vina into PDB 4KD1 (CDK2 monomer, no cyclin), exhaustiveness 8; about 2,000 molecules per policy plus replicate seeds; PoseBusters (22 checks); hinge proximity. Redocking of the native ligand comes first.
* **Plan.** The analysis plan (`campaigns/c2026-10-03/plan/`) was frozen **before this campaign's runs**, but it was written after the project's earlier exploratory work and with those historical results known; it is not a plan written before all exploration. Amendments are appended in `CAMPAIGN_PLAN_AMENDMENTS.md`.

## Counting definitions (read before quoting any number)

* **Population-row statistics** (`13_*` endpoint tables, `05_ga_*` trajectories and figures, `05h` robustness, H5 proximity): rows of the saved population files, **raw stored SMILES, every occurrence**, one value per GA run. This is not the standardized-distinct view the frozen plan names as primary; the endpoint code does not construct it (`KNOWN_ISSUES.md` #2).
* **Docking sets** are drawn as distinct **raw SMILES** within each set; "distinct molecules" in docking figure titles means that.
* `05g_positive_control_summary_v2.csv` saves both a weighted view and a standardized-distinct view (columns suffixed `_weighted`, `_distinct`).
* **Generation 0** in trajectories is the set of starting molecules, which are **training-set molecules** (similarity 1.0; in-sample predictions). Later generations show generated molecules only.
* **"Transfer"** is always one named quantity: surrogate prediction, fingerprint proximity (ECFP4 Tanimoto >= 0.6 to a held-out active), or docking score. Measured activity of generated molecules does not exist.

## What the saved tables show (campaign `c2026-10-03`, per policy; details in `presentation/CLAIMS_TO_EVIDENCE.md`)

Unit of replication is the GA run (n = 5 per arm), so ranges and counts of runs matter more than p-values (the smallest exact two-sided sign-flip p with n = 5 is 0.0625).

* **Surrogate (shared).** Held-out R2 about 0.64 / 0.49 / 0.30 on random / scaffold / paper splits; controls near zero (`shared/results/04_surrogate_metrics.csv`). Scaffold and paper splits are not chemical-series independence or temporal validation.
* **H1 (predicted activity rises).** The median predicted pActivity of generated molecules rises in `multi_real` and `activity_only` in all 5 runs of both policies; every generation step is non-decreasing in all 5 runs only for `activity_only`. This is the surrogate's own prediction.
* **H2 (similarity falls as prediction rises).** For `multi_real`, generation-median similarity to the training set falls as prediction rises in 4 of 5 runs per policy; the drug-likeness-only arm, which has no activity term, ends at a similar similarity, and `activity_only` moves toward training chemistry. Starting similarity is 1.0 for every arm, so ending below it is not specific to an arm. The planned within-size-stratum repetition is not implemented.
* **H3 (docking transfer): not assessed.** The validity gate fixed in the frozen plan (does Vina separate held-out actives from property-matched decoys?) was not passed in any group or policy, so docking-based H3 and the docking-based H5 contingency were not carried out. This is an **inconclusive** benchmark: its intervals include chance and also possibly useful discrimination, so it does not show that Vina lacks signal, and the decoys are putative, not measured, inactives. Docking scores fall with molecular size.
* **H4 (scrambled-label curve).** The scrambled-surrogate arm shows a rising curve of the surrogate it follows (rank correlation >= 0.9 in 5 of 5 runs), while the real surrogate's prediction of its molecules is not above the drug-likeness baseline. "No transfer" here concerns real-surrogate predictions, not biological activity; one fixed label permutation was used.
* **H5 (similarity floor): exploratory and partly provisional.** The paired per-seed change in final predicted activity under a floor is negative in at least 4 of 5 seeds only for `activity_only` at some floors, and not for `multi_real`. The saved exchange-rate estimates and intervals (`13_h5_summary.csv`, `exchange_*`) come from an estimator that fills missing proximity with zero and conditions bootstrap resamples on the denominator's sign; **they are provisional and are not used for any claim**. The saved tables contain exchange intervals of both signs in each policy and I make no statement about them beyond that. Fingerprint proximity shares geometry with the imposed floor and is not an independent endpoint.
* **Held-out proximity.** `activity_only` molecules lie closer in fingerprint space to held-out actives than other arms (`05g_positive_control_summary_v2.csv`). Proximity is neither identity nor activity, and most held-out actives retain a close training neighbour.
* **Docking-pose checks.** Redocking of dinaciclib gives a top-pose RMSD of 0.63-0.65 A over three Vina seeds: pose recovery for this ligand, not an affinity ranking. Most poses pass the 22 PoseBusters checks; that indicates geometric and chemical plausibility, not binding or synthesizability. Hinge proximity (ligand N/O within 3.5 A of the hinge backbone; distance only) differs between arms and is exploratory.
* **Preparation policy.** The two policies give qualitatively similar results in the tables above. They differ in validity: radical-bearing candidates appear only under `legacy` (and `no_radicals` PoseBusters failures only there). The cause cannot be assigned to one of the bundled changes.

## Limitations

Assay heterogeneity (ATP concentration mostly unrecorded; many assays on CDK2/cyclin A while docking uses monomeric CDK2); censored `>` records dropped and records aggregated by median (modelling choices with unquantified selection effects); no measurable noise floor;
ECFP4 cannot separate molecules with identical fingerprints (e.g. stereoisomers) and similarity-based models tend to smooth over activity cliffs; the GA's edits are chemically naive and QED rewards resemblance to existing drugs; Vina is a rigid-receptor score; decoys are putative inactives and some descriptors remain imbalanced;
"held-out" actives are neither distant chemistry nor, in general, free of training siblings; the scrambled-surrogate runs vary the search for one label permutation, and the size-only control is a regression, not a GA control; n = 5 runs per arm; nothing was validated experimentally.
Redocking, docking discrimination, geometric plausibility, hinge proximity, synthetic feasibility and biological validation are different things and only the first and the plausibility checks were available here. An inconclusive comparison is not evidence of no effect, and failing to show superiority is not equivalence.

## Repository map

| Where | What |
|---|---|
| `campaigns/c2026-10-03/` | The current campaign: `inputs/`, `plan/` (frozen plan + amendments), `manifests/` (ledger, hashes, logs), `shared/`, `legacy/`, `corrected/` (results and figures per policy), `docking_store/`, `superseded/`, `HANDOFF.md` (written before the review; superseded where it differs from this README and `KNOWN_ISSUES.md`) |
| `presentation/` | Selected figure copies, captions, claims table, start page |
| `archive/` | 123 historical (pre-campaign) files moved from `data/`, `results/`, `figures/` with checksums, plus pre-cleanup copies of documents; `RESTORE.md` |
| `historical_record/` | Pre-campaign hash manifests and copies of the earlier README and CLAUDE.md |
| `src/cdk2moo/`, `scripts/`, `tests/` | Code (modules hold logic; scripts are numbered stages that read and write files) |
| `data/raw/`, `data/structures/` | The three frozen campaign inputs at their original paths (ChEMBL CSV, 4KD1, crystal ligand) |
| `environment.yml`, `pyproject.toml`, `reproducibility/` | Dependency specification, packaging, environment snapshot |
| `KNOWN_ISSUES.md`, `decisions.md`, `CLAUDE.md` | Open problems; decision log (D-45 describes this cleanup); project rules |

## Reproduction (instructions for a future run; not executed or verified by this cleanup)

From the repo root, in the `cdk2` conda environment (`environment.yml`; the package is expected to be installed in editable mode, and the tests need the scientific packages, not only the `unittest` runner):

```
python -m unittest discover -s tests
python scripts/campaign_run.py init   --campaign <new-campaign-id>
python scripts/campaign_run.py run    --campaign <new-campaign-id> --cores 8
python scripts/campaign_run.py status --campaign <new-campaign-id>
```

A new campaign id repeats the experiment with its own plan; do not re-initialise `c2026-10-03` over existing outputs. The recorded campaign took roughly 8.5 h of docking on 8 cores. `campaign_run.py verify-historical` expects the historical files at their original paths: restore them first (`archive/RESTORE.md`).
Test results are attributed to their recorded date: 159 tests passed on 2026-10-05 in the review's run (`../review-2026-10-05/tests.log`); they were not rerun in this cleanup. Resume checking is weaker than documented (`KNOWN_ISSUES.md` #5), and the recorded Git commit does not contain the reviewed source (#8), so cloning the repository is not sufficient to reproduce this project.

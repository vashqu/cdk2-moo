# Hand-off for audit: campaign c2026-10-03 (2026-10-05)

Everything below is traceable to a file under `campaigns/c2026-10-03/` unless it names another path. Historical (pre-campaign) artifacts are untouched; nothing here pools policies, groups or experiments.

## 1. Implementation changes (this task; historical outputs unchanged)

| Area | Files | What changed |
|---|---|---|
| Campaign isolation | `src/cdk2moo/config.py`, `campaign.py`, `campaign_driver.py`, `scripts/campaign_run.py`, `scripts/campaign_inventory.py` | Campaign/scope paths from environment variables; no campaign -> stop; historical only via `CDK2_CAMPAIGN=historical`; hashing, ledger, stage manifests with consumed-upstream hashes, source/environment identity, 69-stage dependency table and scheduler |
| Docking integrity | `dock_validate.py`, `dock_jobs.py`, `docking.py`, `scripts/08d_dock_sets.py` | Atomic per-job pose + record, SHA-256 of pose, job-id chain, pose identity (InChIKey connectivity; stereo when fully specified), schema, writer lock, quarantine instead of deletion, scoped repair. In-campaign fix: raw-string `job_id` (A5) |
| Distinct-molecule metrics | `recovery.py`, `docking_stats.py`, `endpoints.py`, `hinge.py`, `scripts/13_hypothesis_endpoints.py` | Order-invariant distinct-molecule metrics primary, weighted secondary; cluster bootstrap; H1-H5 endpoints as frozen |
| Policies and matched starts | `candidates.py`, `starts.py`, `ga.py`, `mutations.py`, `standardize.py`, `scripts/03b_eligible_starts.py`, `05_run_ga.py`, `05f_*`, `10_run_constrained_ga.py` | Legacy vs corrected preparation; one eligible-starts list; similarity floor as hard feasibility filter; commit_csv/manifest protected outputs (`run_outputs.py`) |
| Docking sets and analysis | `docksets.py`, `scripts/08c`, `08e`, `08f`, `08g`, `12_*` | Groups P/M/S with stated comparator role; validity gate; PoseBusters 22 checks, three states; score-to-pose coverage |
| Figures/tables | `figures_model.py`, `figures_validation.py`, `scripts/11_make_figures.py`, `11b_table1.py` | Titles state campaign, policy, n, weighting, uncertainty; fig5 per group; fig6 pose checks; group-aware Table 1 |
| Code style | scripts | No `def` other than `main` in `scripts/` (helpers moved to `src/`) |

## 2. Tests

`python -m unittest discover -s tests` -> **Ran 159 tests, OK** (2026-10-05, after the last code change). Tests establish that metrics, stores, validators and the campaign mechanism behave as specified (including a
historical-equivalence test using explicit historical paths and a 13-test campaign integration test that checks no stage falls back to historical files). They do not establish that an experiment answers its scientific question.

## 3. Frozen plan

`plan/CAMPAIGN_PLAN.md` (sha256 prefix 18e4d79245021d49, recorded in `plan/plan.sha256.json` and as the ledger event `plan_frozen`), byte-identical today. Amendments are appended only in `plan/CAMPAIGN_PLAN_AMENDMENTS.md` (A1-A7).
The plan was written knowing the historical results; items shaped by that knowledge are marked in the plan.

## 4. Execution ledger (`manifests/ledger.jsonl`)

70 stage entries: 69 planned stages **completed** + `plan_frozen`; 82 `running` events, 71 `completed` events (69 stages; the two `figures` stages completed twice because the figures were regenerated after the title fix, A6), **11 failed attempts**, all repaired and rerun (none removed):
`holdout_eval` x4 (A4), `legacy/dock` x1 (A5), `table1` x2 and `figures` x4 (A6). Each completed stage has a manifest (`manifests/stages/*.json`) with output hashes, consumed-upstream hashes, log hash and the source-tree hash it ran under.
Stage logs: `manifests/logs/`. Failed outputs set aside: `superseded/`. Docking store: 3,660 job records (3,654 ok, 6 prep_failed, all recorded), 0 anomalies on re-validation.

## 5. Manifests and preservation

* Inputs: `manifests/inputs.json` (cached ChEMBL download e42b76cf..., 4KD1 f902a651..., crystal ligand 154ef0db...; each equals the historical hash).
* Source: `manifests/source_identity_at_init.json` (c7f4d920...) vs `manifests/final_inventory.json` (9d161d7b...; 7,768 files hashed); differences listed in A7. Environment: `manifests/environment_at_init.json` (new capture; says nothing about the environment of the historical results).
* Historical preservation: `historical_record/historical_hashes_2026-10-03_pre-campaign.json`; `python scripts/campaign_run.py verify-historical` -> **132 recorded files; changed 0, missing 0, added since 0**. The only historical-record files that differ now are the three the task required rewriting: `decisions.md` (append-only: its first 89,380 bytes still hash to the recorded value; entries D-40..D-44 appended), `README.md` and `CLAUDE.md`. Copies whose hashes equal the recorded pre-campaign hashes: `historical_record/README_before_campaign_c2026-10-03.md` and `historical_record/CLAUDE_before_campaign_c2026-10-03.md` (the latter reconstructed by reversing the edits and checked against the recorded hash). `verify-historical` therefore reports `changed 3` (those three), `missing 0`.

## 6. New figures and tables

Per policy: `<policy>/figures/final/` fig1_surrogate, fig2_trajectories, fig3_pareto_fronts, fig3b_pareto_structures, fig4_similarity_vs_prediction, fig5_docking_{P,M,S}, fig6_pose_checks, table1.md; `<policy>/results/table1.csv`.
Key tables: `results/08e_validity_gate.csv`, `08e_h3_endpoints.csv`, `13_h1_summary.csv`, `13_h2_summary.csv`, `13_h4_summary.csv`, `13_h5_summary.csv`, `05g_positive_control_summary_v2.csv`, `05h_robustness_summary.csv`, `08g_hinge_by_set.csv`, `audit/`.

## 7. Claim-to-evidence table (legacy / corrected values; n = 5 GA runs unless stated)

### Campaign-supported findings

| Claim | Experiment / policy | Source artifact | Metric and denominator | Uncertainty | Controls | Limitations | README |
|---|---|---|---|---|---|---|---|
| Median real prediction rises in activity-driven arms (H1 "rises") | Primary GA, both policies | `results/13_h1_summary.csv`, `13_h1_runs.csv` | net rise m50-m0: multi_real +1.07 / +1.13; activity_only +1.91 / +1.93; 5/5 runs | min-max over runs (1.02-1.21; 1.03-1.40) | druglike_only -0.51 / -0.52; scrambled -0.53 / -0.65 (median real prediction) | The surrogate's own prediction | Results table, H1 |
| "Monotonically" holds only for activity_only | same | same | runs with every step non-decreasing: 5/5 activity_only; 0/5 multi_real (86% of steps) | counts of runs | - | Definition fixed in plan | H1 |
| H2 literal holds for multi_real but is not attributable to the activity term | Primary GA | `13_h2_summary.csv` | Spearman(generation, median similarity): median -0.95 / -0.92, negative in 4/5; multi_real minus druglike_only final similarity +0.024 / +0.007 (farther in 1/5, 2/5 seeds) | min-max over runs | druglike_only drifts to 0.36 | Raw ECFP4 similarity depends on molecular size (size-conditioned versions in the same files) | H2 row |
| activity_only moves toward training chemistry, not away | Primary GA | `13_h2_summary.csv`; `05h_robustness_summary.csv` | final median similarity ~0.70, +0.35 vs druglike_only; negative Spearman in 0/5 runs; stable in 5/5 scaffold splits | range over runs and splits | druglike_only baseline | Large, low-QED molecules (median 35 atoms, QED ~0.3) | H2 row; descriptive paragraph |
| H4 curve similar, no transfer | Primary GA | `13_h4_summary.csv` | Spearman >= 0.9 in 5/5 runs both arms; real-surrogate prediction of multi_scrambled minus druglike_only -0.02 / -0.13 vs multi_real +1.56 / +1.64 | min-max over runs; scrambled above baseline in 2/5 seeds | scrambled labels (fixed permutation seed 20260930), druglike_only | One scrambled permutation | H4 row |
| Vina does not separate held-out actives from matched decoys (gate failed) | Docking, groups M and S (held-out), P (in-sample), both policies | `results/08e_validity_gate.csv` | AUC(-Vina; actives vs decoys), n actives/decoys 100/100 (M), 93/93 (S); lower bound < 0.5 raw and stratified everywhere | 95% cluster-bootstrap, 2,000 resamples, clusters = GA seed (decoys) / scaffold (actives) | matched decoys, random molecules, size strata, 3-seed engine noise (median SD 0.012 kcal/mol) | Rigid receptor; monomeric CDK2; activity data largely CDK2/cyclin A | H3 row |
| Docking score tracks heavy-atom count | All docked molecules per group | `08e_docking_analysis` log, `fig5_docking_*.png` | -0.04 to -0.07 kcal/mol per heavy atom; Spearman -0.3 to -0.5 | descriptive | - | Linear fit over all molecules of the group | H3 row |
| Corrected policy removes radical-bearing candidates (legacy has them) | Primary GA; docking audit | `results/05_ga_populations.csv` both scopes; `audit/` | 149 of 12,095 distinct (arm, SMILES) legacy generated molecules (activity_only 106) vs 0 / 12,030; 8 legacy poses fail `no_radicals` | counts | policy comparison | Combined intervention, not decomposed | Preparation paragraph |
| Same verdicts under both policies | H1-H5 rows | files above | see tables | - | - | n = 5 per arm | Results |
| H5 attenuation exists only for activity_only at some floors; multi_real constrained runs lose no predicted activity | Constrained GA, 3 experiments | `13_h5_summary.csv` | paired final median difference, constrained minus twin; attenuation if negative in >= 4/5 seeds | min-max over seeds | unconstrained twins matched by seed and start | Floors are one operationalization of "stay near training" | H5 row |
| Surrogate accuracy and controls | Surrogate evaluation, shared | `shared/results/04_surrogate_metrics.csv` | median over 10 split seeds: R2 0.64/0.49/0.30 (random/scaffold/paper) vs scrambled -0.12/-0.13/-0.08 | min-max in file | scrambled, mean, size-only | Overlapping test sets (descriptive) | Surrogate paragraph |
| Redocking passes | Receptor/redock, shared | `shared/results/08b_redocking.json` | top-pose RMSD 0.63-0.65 A, 3 seeds (< 2 A) | 3 seeds | crystal ligand | Says nothing about activity | Docking checks |

### Historical-only (not re-established in this campaign)

Historical docking, H5 and positive-control numbers in the pre-campaign README (`historical_record/README_before_campaign_c2026-10-03.md`, `results/`, `figures/final/`) are historical; where this campaign re-measured the same quantity the campaign value is used above. The historical "3 of 6 claims stable" became 4 of 6 under both policies here.

### Exploratory (not used for any hypothesis statement)

Hinge proximity by arm (activity_only 76-92% vs druglike_only 38-55%, held-out actives 88-95%; distance only, not a hydrogen bond); multi_real's collapse to the window minimum (median 21 heavy atoms) and activity_only's growth (35); fingerprint proximity of generated molecules to held-out actives (activity_only 84% / 90% at >= 0.6, molecule level; 88% of held-out actives already have a training sibling at >= 0.6); P_sup values in groups M, S, P (tabulated, descriptive, not interpretable as transfer); exchange rates (defined with an interval excluding 0 only for scaffold activity_only at floor 0.7 in both policies).

### Unanswered

Whether any generated molecule is active; H3 and the H5 transfer claim (no validated orthogonal signal); the contribution of each of the three components of the corrected policy; whether "held-out" chemistry is unseen (siblings remain in training); temporal generalization; a random-ChEMBL reference.

## 8. Reproduction commands

See README "Reproduce". Minimal: `python scripts/campaign_run.py init --campaign c2026-10-03 && python scripts/campaign_run.py run --campaign c2026-10-03 --cores 8`, then `python -m unittest discover -s tests` and
`python scripts/campaign_run.py verify-historical`. Expect the docking stages to take about 8.5 h on 8 cores. A new campaign id repeats the experiment with its own plan; `c2026-10-03` must not be re-initialised over existing outputs.

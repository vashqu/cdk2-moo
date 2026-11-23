# Hand-off: campaign c2026-11-22 (2026-11-23)

Status: **complete**, 38 of 38 stages, no failed attempts in the ledger (`manifests/ledger.jsonl`). Ran on commit `a48134a`, plan hash `9aa550bf8b31…` (`plan/plan.sha256.json`), started 2026-11-22 21:37 UTC, last stage 2026-11-23 11:22 UTC.

## What was run
`corrected` policy only; 10 GA seeds (42-51); one scrambled-label permutation per seed (permutation seed = 20260930 + seed - 42); distinct-standardized counting enforced by `cdk2moo.distinct`; element-wise stereo check on every docked pose; repaired H5 estimator; resume that re-checks consumed inputs and the code hash.

## Results
See `README.md` (results table) and `docs/claims_to_evidence.md` (rows A1-A7). Tables: `corrected/results/13_*.csv`, `08e_*.csv`, `table1.csv`. Figures: `corrected/figures/final/`.

## Deviations and failures
1. Within-heavy-atom-stratum repetition of H2 was required by the plan and not produced (`plan/CAMPAIGN_PLAN_AMENDMENTS.md` P1; `docs/known_issues.md` #12).
2. 4 poses rejected by the stereo check and 1 preparation failure out of 2,200 docking jobs (`corrected/results/08e_denominators.csv`).
3. A kill-and-resume during docking was not tested.

## Not in git
`corrected/results/10_constrained_{primary,molecule,scaffold}_populations.csv` (about 120 MB each) exceed GitHub's 100 MB limit and are git-ignored; they exist locally and their hashes are in `manifests/stages/corrected__constrained_*.json`.

## Not done
No LICENSE file. The campaign folder is about 975 MB (the docking store dominates); decide before publishing whether to track it in git. No claim here is about activity of any generated molecule.

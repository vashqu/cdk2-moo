# Amendments to the frozen plan of campaign c2026-10-03

`CAMPAIGN_PLAN.md` is byte-identical to the version hashed at freeze (`plan.sha256.json`). Everything below is a dated, append-only note. None of it changes a
seed, a parameter, a threshold or a decision rule in the plan, and none of it was prompted by a campaign outcome: items A1-A3 were made before any analysis output was
read; A4-A6 are repairs to code that failed or produced a defective artifact, made without any change to what is computed.

**A1 (2026-10-04, before any campaign result was read). Stage count.** The plan says 67 stages. The driver has 69: the endpoint analysis (`scripts/13_hypothesis_endpoints.py`,
stage `endpoints`, per policy) was written after the freeze to compute the H1-H5 quantities exactly as defined in section 3; it adds no new endpoint.

**A2 (before any result). Docking group tables.** `08e_docking_analysis.py` and `11b_table1.py` were adapted to the group (P, M, S) structure of section 2; `11_make_figures.py` was
written for campaign scopes; titles of figures state campaign, policy scope, sample sizes, weighting and uncertainty.

**A3 (before any result). Docking-set rules.** `08c_select_docking_sets.py` rules (100 per set, keyed seeds, matched decoys) were fixed before it ran on campaign data; a nested helper
was inlined in it (no change of behaviour).

**A4 (after results existed). `05g_positive_control_eval.py` crashed twice** (stage `holdout_eval`, both policies): it read `median_proximity_to_held_actives` without the `_weighted`
suffix and unpacked two values from a function that returns three. Fixed; the stage was rerun in full. The failed attempts are in the ledger.

**A5. Pose validator false positive.** `dock_validate.check_pose` read the pose's `job_id` through RDKit `GetPropsAsDict`, which converts an id made only of digits (or digits with one
"e") into a number, so two correct poses (jobs 4445476191389595 and 84868095015190e9) were refused (`pose_mismatch`) and `legacy/dock` failed at export. The check now reads raw strings;
both poses have the same InChIKey as their jobs; a regression test was added. No score, pose or job was changed; the store re-validated with no anomalies.

**A6. Output declarations and titles.** Stage `figures` declared output names that the script does not write, and the fig 1/2/4 titles ran off the figure; declarations were corrected and
titles wrapped; the first figure set is kept in `superseded/`. `table1` stage: `11b_table1.py` was rewritten for group-aware sets (A2).

**A7. Code identity.** Source hash at init: `c7f4d920...`; at the final inventory: `9d161d7b...`. Files changed since init: scripts 05g, 07c, 08c, 08e, 08f, 08g, 11, 11b, 12; src
alerts.py, campaign_driver.py, dock_validate.py, figures_model.py, figures_validation.py; tests/test_dock_jobs.py; added: scripts/13, src docking_stats.py, endpoints.py, hinge.py and three
test files. No file used by curation, features, splits, surrogates, GA, constrained GA or receptor preparation changed. The source hash that each stage ran under is in its manifest
(`manifests/stages/*.json`, field `source`).

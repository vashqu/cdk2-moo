# Amendments to the plan of campaign c2026-11-22

(none yet)

## Post-run record (2026-11-23): required output not delivered

**P1. Within-heavy-atom-stratum repetition of H2 was not produced.** Section 3 of the frozen plan makes it a required output of stage 13 and says its absence is "reported as a failure of the stage". `scripts/13_hypothesis_endpoints.py` computes the size-conditioned percentile version (`delta_size_percentile_median`) and the heavy-atom difference (`delta_heavy_atoms_median`), but has no per-stratum computation. Stage 13 is recorded as completed; this is the stated failure. The H2 "attributable" result below is therefore reported with the size-percentile repetition only. Implementing the missing computation changes `scripts/`, which re-runs every stage under the code-hash rule; a post-hoc analysis would have to be labelled as such.

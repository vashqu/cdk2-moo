# Start here

This folder is an **organized presentation package, not a newly validated release**. The full scientific record is the rest of the project and `archive/`. Nothing here was recomputed in the cleanup of 2026-10-05, and clean-checkout reproducibility has not been verified.

1. **The question and methods, current results, limitations:** [`../README.md`](../README.md) (campaign `c2026-10-03`, two preparation policies, never pooled).
2. **Evidence:** [`CLAIMS_TO_EVIDENCE.md`](CLAIMS_TO_EVIDENCE.md) (claim, policy, source table/columns, counting definition, qualification) and [`FIGURE_INDEX.md`](FIGURE_INDEX.md) (captions, provenance, limitations; copies in `figures/`).
3. **Unresolved issues:** [`../KNOWN_ISSUES.md`](../KNOWN_ISSUES.md). In short: the H5 exchange-rate estimator is unreliable (its outputs are provisional); the endpoint pipeline counts raw population rows, not standardized distinct molecules as the plan intended; stereochemistry validation of docked poses has gaps; constrained-run provenance was discarded; resume validation is incomplete; the committed Git history does not contain the reviewed source. These were **not repaired** in this pass.
4. **History and decisions:** [`../decisions.md`](../decisions.md) (append-only; D-45 describes this cleanup), [`../archive/ARCHIVE_INDEX.md`](../archive/ARCHIVE_INDEX.md), [`../CLAUDE.md`](../CLAUDE.md) (project rules; original hypotheses kept as history).
5. **Independent review:** `../../review-2026-10-05/REVIEW.md`.

What the project does and does not show (details and qualifications in the README): it shows how optimization curves behave under a learned objective, a scrambled control and a drug-likeness-only baseline; it did **not** establish whether any generated compound is active, did **not** measure the surrogate's error on generated molecules, and its docking benchmark did not establish discrimination (inconclusive, not shown to be absent).

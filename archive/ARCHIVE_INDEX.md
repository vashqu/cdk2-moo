# Archive index

Created 2026-10-05 in a documentation-and-archiving pass. Nothing was recomputed, regenerated or deleted. File operations were limited to copying, moving and checksumming.

| Directory | Contents | Status |
|---|---|---|
| `historical_pre_campaign/` | 123 files from the top-level `data/`, `results/`, `figures/` (the pre-2026-10-03 project: historical populations, surrogate and docking tables, poses, receptor intermediates, figures), original relative paths preserved | Moved; checksums before = after for all 123 (`archive_manifest.csv`) |
| `documentation_before_cleanup/` | Copies of `README.md`, `CLAUDE.md`, `decisions.md`, `campaigns/c2026-10-03/HANDOFF.md` as they stood before this pass | Copies; the originals stay in place. `decisions.md` is append-only and was only appended to |
| `archive_manifest.csv` | One row per archived or copied file: original path, archive path, bytes, SHA-256 before and after, status | 123 `move` + 4 `copy_document` rows, all verified |
| `RESTORE.md`, `restore_historical.py` | How to put files back without overwriting | Script written, **not run** |

## What is historical, and why it was archived

The files are the output of the project **before** campaign `c2026-10-03`: an earlier set of GA runs, surrogate evaluation, docking (788-molecule panel) and figures, including the relaxed-window pilot and
definitions later corrected (`_v2` files). They are historical evidence, not the current evidence: counts, definitions and verdicts differ from the campaign, so they are kept out of the working tree to avoid mixing the two. The review of 2026-10-05 recommended archiving them intact
(`review-2026-10-05/FILE_RETENTION.csv`, recommendation `ARCHIVE_HISTORICAL`); this pass checked that every listed file existed on disk and that the repository had not changed since the review before moving anything.
Their original hashes are still recorded in `historical_record/` (not rewritten). One archived file, `data/raw/.gitkeep`, is an empty placeholder.

## What was deliberately NOT moved

* `data/raw/chembl_CHEMBL301_activities.csv`, `data/structures/4KD1.pdb`, `data/structures/4kd1_1QK_crystal.sdf`: `init_campaign` reads these original paths.
* Everything under `campaigns/c2026-10-03/` (both policy scopes, constrained outputs, docking store, maps, exports, logs, superseded attempts, inventories). The review found all 7,770 inventoried files matched their hashes; this pass did not touch them.
* `historical_record/`, `src/`, `scripts/`, `tests/`, `environment.yml`, `pyproject.toml`, `setup.py`, `reproducibility/`, `notebooks/`, `.git/`.
* `__pycache__/` directories (recommended for removal by the review; left alone because cache cleanup was optional and nothing needs it removed).

## Consequences

`verify-historical` will report the 123 archived files as missing until they are restored; see `RESTORE.md`. `results/`, `figures/` and `data/` each contain a `RELOCATED.md` pointing here.

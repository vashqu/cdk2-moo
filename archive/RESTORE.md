# Restoring archived historical files

What was moved: 123 files from the top-level `data/`, `results/` and `figures/` into `archive/historical_pre_campaign/`, **keeping their original relative paths** (for example
`results/05_ga_populations.csv` is now `archive/historical_pre_campaign/results/05_ga_populations.csv`). Every move is in `archive_manifest.csv` with its size and SHA-256 before and after the move (all 123 identical).
Not moved: `data/raw/chembl_CHEMBL301_activities.csv`, `data/structures/4KD1.pdb`, `data/structures/4kd1_1QK_crystal.sdf` (campaign initialisation reads these paths), everything in `campaigns/`, `historical_record/`, source, tests.

## Restore (no overwriting)

```
python3 archive/restore_historical.py            # dry run: prints what would be copied and any conflicts
python3 archive/restore_historical.py --apply    # copies archive files back to their original paths
```

The script copies (the archive stays complete), checks each archive file against the recorded checksum first, leaves an existing identical file alone, and skips and reports any existing file with different content. It was written during the cleanup and has **not been run**.
Manual equivalent for one file: `cp -n archive/historical_pre_campaign/<path> <path>` (`-n` refuses to overwrite).

## What changes when files are archived

* `python scripts/campaign_run.py verify-historical` compares against the original paths in `historical_record/`; with the files archived it will report them as **missing**, not as altered. Restore, then verify, to check them. (Before the cleanup it reported 0 missing and 3 changed: README.md, CLAUDE.md and the append-only decisions.md.)
* `CDK2_CAMPAIGN=historical` analysis and `tests/test_legacy_equivalence.py` (which skips itself when the historical populations are absent) need the files restored.
* The current campaign `campaigns/c2026-10-03/` is unaffected; it does not read the archived files. Nothing inside it was moved or edited.
* The hash manifests in `historical_record/` were not rewritten.

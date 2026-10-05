#!/usr/bin/env python3
"""
Put archived historical files back at their ORIGINAL paths (so `campaign_run.py verify-historical` and CDK2_CAMPAIGN=historical can find them).

Reads archive/archive_manifest.csv. For every `move` row it checks the archive copy against the recorded checksum, then copies it back. It never overwrites:
if the original path already holds a file with the same checksum it is left alone, and if it holds DIFFERENT content the file is skipped and reported.
Default is a dry run; pass --apply to copy. Copies, not moves: the archive stays complete.

Run from the repo root:   python3 archive/restore_historical.py            (dry run)
                          python3 archive/restore_historical.py --apply
"""

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--apply", action="store_true")
args = ap.parse_args()
root = Path(__file__).resolve().parent.parent
counts = {"restore": 0, "already_present": 0, "conflict_skipped": 0, "archive_checksum_bad": 0}
for row in csv.DictReader(open(root / "archive" / "archive_manifest.csv")):
    if row["action"] != "move":
        continue
    source, target = root / row["archive_path"], root / row["original_path"]
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != row["sha256_after"]:
        print("ARCHIVE CHECKSUM DIFFERS:", row["archive_path"])
        counts["archive_checksum_bad"] += 1
        continue
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest() == digest:
            counts["already_present"] += 1
        else:
            print("CONFLICT, not overwritten:", row["original_path"])
            counts["conflict_skipped"] += 1
        continue
    counts["restore"] += 1
    if args.apply:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
print(counts, "(applied)" if args.apply else "(dry run: nothing copied)")

#!/usr/bin/env python3
"""
Record (and later verify) the hashes of all historical artifacts, before any campaign work edits or adds anything.

Writes : historical_record/historical_hashes_<label>.json   every file in data/, results/, figures/, reproducibility/,
                                                              decisions.md, CLAUDE.md, README.md, environment.yml
         historical_record/source_identity_<label>.json     content hash of src/, scripts/, tests/ and the Git status
Verify : python scripts/campaign_inventory.py --verify historical_record/historical_hashes_<label>.json

Files are tagged "historical" (older than the earlier repair pass) or "repair_pass" (new files that pass wrote); both are
preserved. This script does not need a campaign selection: it only reads and hashes.
"""

import argparse
import json
from pathlib import Path

from cdk2moo import campaign, config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="2026-10-03_pre-campaign")
    ap.add_argument("--verify", default=None, help="path of a recorded inventory to compare the current files with")
    args = ap.parse_args()
    root = config.PROJECT_ROOT

    if args.verify:
        recorded = json.loads(Path(args.verify).read_text())
        result = campaign.verify_inventory(root, recorded["files"])
        print(f"{len(recorded['files'])} recorded files; changed {len(result['changed'])}, missing {len(result['missing'])}, "
              f"added since {len(result['added'])}")
        for key in ("changed", "missing"):
            for rel in result[key][:20]:
                print(f"  {key}: {rel}")
        raise SystemExit(1 if result["changed"] or result["missing"] else 0)

    out = root / "historical_record"
    out.mkdir(exist_ok=True)
    inventory = campaign.historical_inventory(root)
    target = out / f"historical_hashes_{args.label}.json"
    if target.exists():
        raise SystemExit(f"{target} exists; historical records are never overwritten")
    tags = {}
    for info in inventory.values():
        tags[info["tag"]] = tags.get(info["tag"], 0) + 1
    campaign.write_json_atomic(target, {"recorded_utc": campaign.now_utc(), "label": args.label, "tag_counts": tags,
                                        "note": "historical = modified before 2026-10-03 00:00; repair_pass = written by the "
                                                "2026-10-03 repair pass. No manifests existed for any of these files.",
                                        "files": inventory})
    campaign.write_json_atomic(out / f"source_identity_{args.label}.json", campaign.source_identity(root))
    print(f"hashed {len(inventory)} files {tags}; wrote {target.name}")


if __name__ == "__main__":
    main()

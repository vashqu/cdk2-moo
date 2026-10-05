#!/usr/bin/env python3
"""
Initialise, run, and inspect a campaign: the full repaired pipeline in campaigns/<id>/, kept separate from the historical results.

  python scripts/campaign_run.py init   --campaign <id>                 create the tree, freeze and hash the inputs, plan every stage
  python scripts/campaign_run.py run    --campaign <id> [--only 'shared/*' 'legacy/ga_*'] [--cores 8]
  python scripts/campaign_run.py status --campaign <id>                 latest ledger status of every stage
  python scripts/campaign_run.py verify-historical                      check the historical artifacts against the recorded hashes

A stage runs only when every dependency has a complete manifest with unchanged output hashes. Re-running `run` resumes: completed stages are
skipped, failed or interrupted ones are redone (their old outputs are moved to superseded/, never deleted). The stage table, the ledger and the
manifests are described in cdk2moo/campaign_driver.py. Nothing here needs a CDK2_CAMPAIGN environment variable; the driver sets it for each stage.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from cdk2moo import campaign, campaign_driver as driver
from cdk2moo.config import PROJECT_ROOT

HISTORICAL = PROJECT_ROOT / "historical_record" / "historical_hashes_2026-10-03_pre-campaign.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["init", "run", "status", "verify-historical"])
    ap.add_argument("--campaign", default=None)
    ap.add_argument("--only", nargs="+", default=None, help="stage patterns such as 'shared/*' or 'legacy/ga_*'")
    ap.add_argument("--cores", type=int, default=8)
    args = ap.parse_args()

    if args.command == "verify-historical":
        sys.exit(subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "campaign_inventory.py"), "--verify", str(HISTORICAL)]).returncode)
    if not args.campaign:
        sys.exit("--campaign <id> is required")
    if args.command == "init":
        inputs = driver.init_campaign(PROJECT_ROOT, args.campaign, HISTORICAL)
        print(f"campaign {args.campaign} initialised; {len(driver.STAGES)} stages planned")
        for rel, info in inputs.items():
            print(f"  frozen input {rel}: {info['sha256'][:16]}  matches historical record: {info['matches_historical_record']}")
    elif args.command == "status":
        latest, counts = driver.status(PROJECT_ROOT, args.campaign)
        print(json.dumps(counts))
        for key in driver.BY_KEY:
            print(f"  {latest.get(key, 'unplanned'):<30} {key}")
    else:
        failed = driver.run_campaign(PROJECT_ROOT, args.campaign, args.cores, args.only)
        latest, counts = driver.status(PROJECT_ROOT, args.campaign)
        print(f"stage status counts: {counts}")
        if failed:
            sys.exit(f"FAILED stages: {failed}; see campaigns/{args.campaign}/manifests/logs/")


if __name__ == "__main__":
    main()

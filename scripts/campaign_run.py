#!/usr/bin/env python3
"""
Initialise, run, and inspect a campaign: the full repaired pipeline in campaigns/<id>/, kept separate from the historical results.

  python scripts/campaign_run.py init   --campaign <id>                 create the tree, freeze and hash the inputs, plan every stage
  python scripts/campaign_run.py run    --campaign <id> [--only 'shared/*' 'legacy/ga_*'] [--cores 8]
  python scripts/campaign_run.py status --campaign <id>                 latest ledger status of every stage
  python scripts/campaign_run.py freeze-plan --campaign <id>            hash campaigns/<id>/plan/CAMPAIGN_PLAN.md; `run` refuses without it
  python scripts/campaign_run.py verify-historical                      check the historical artifacts against the recorded hashes

`run` also refuses to start unless (1) the plan is frozen and unchanged and (2) the code in src/ and scripts/ is committed, so the recorded git commit
identifies the code that produced the results (`--allow-uncommitted` is for development runs and says so in the ledger).

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

HISTORICAL = PROJECT_ROOT / "archive" / "historical_record" / "historical_hashes_2026-10-03_pre-campaign.json"


CODE_PATHS = ["src", "scripts", "tests", "pyproject.toml", "environment.yml"]


def uncommitted_code():
    """Lines of `git status` for the code and environment files (untracked files included); empty if everything is committed."""
    out = subprocess.run(["git", "status", "--porcelain", "--", *CODE_PATHS], cwd=PROJECT_ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        return [f"git status failed: {out.stderr.strip()}"]
    return out.stdout.splitlines()


def plan_paths(cid):
    plan_dir = Path(PROJECT_ROOT) / "campaigns" / cid / "plan"
    return plan_dir / "CAMPAIGN_PLAN.md", plan_dir / "plan.sha256.json"


def freeze_plan(cid):
    plan, record = plan_paths(cid)
    if not plan.is_file():
        sys.exit(f"{plan} does not exist; write the plan before freezing it")
    if record.exists():
        sys.exit(f"{record} exists; a frozen plan is never re-frozen (record changes as amendments in CAMPAIGN_PLAN_AMENDMENTS.md)")
    campaign.write_json_atomic(record, {"file": plan.name, "frozen_utc": campaign.now_utc(), "sha256": campaign.sha256_file(plan)})
    print(f"plan frozen: {campaign.sha256_file(plan)[:16]}")


def plan_is_frozen(cid):
    plan, record = plan_paths(cid)
    if not record.exists() or not plan.is_file():
        return False, "no frozen plan (write campaigns/<id>/plan/CAMPAIGN_PLAN.md, then run freeze-plan)"
    if json.loads(record.read_text())["sha256"] != campaign.sha256_file(plan):
        return False, "CAMPAIGN_PLAN.md differs from its frozen hash; record changes in CAMPAIGN_PLAN_AMENDMENTS.md instead"
    return True, ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("command", choices=["init", "run", "status", "freeze-plan", "verify-historical"])
    ap.add_argument("--allow-uncommitted", action="store_true", help="development only: run with uncommitted code (recorded in the ledger)")
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
    elif args.command == "freeze-plan":
        freeze_plan(args.campaign)
    else:
        frozen, why = plan_is_frozen(args.campaign)
        if not frozen:
            sys.exit(f"refusing to run: {why}")
        dirty = uncommitted_code()
        if dirty and not args.allow_uncommitted:
            sys.exit("refusing to run: code is not committed, so the recorded git commit would not identify it:\n  " + "\n  ".join(dirty[:15]))
        if dirty:
            campaign.ledger_append(driver.ledger_path(PROJECT_ROOT, args.campaign), stage="campaign", status="run_with_uncommitted_code", lines=dirty)
        failed = driver.run_campaign(PROJECT_ROOT, args.campaign, args.cores, args.only)
        latest, counts = driver.status(PROJECT_ROOT, args.campaign)
        print(f"stage status counts: {counts}")
        if failed:
            sys.exit(f"FAILED stages: {failed}; see campaigns/{args.campaign}/manifests/logs/")


if __name__ == "__main__":
    main()

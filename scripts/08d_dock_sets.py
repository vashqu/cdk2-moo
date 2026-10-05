#!/usr/bin/env python3
"""
Stage 8d: dock every selected molecule with Vina, resumably and safely, into a docking store SHARED by the campaign's scopes.

Reads : results/08c_docking_sets.csv (this scope's sets), structures/receptor.pdbqt and receptor_box.txt
Writes: <store>/ (default: the campaign's docking_store/)   manifest.json, maps/, jobs/<id>.json|.sdf, quarantine/   (see cdk2moo/dock_jobs.py)
        results/08d_export/docking_scores.csv, poses_seed<N>.sdf      this scope's molecules only (legacy-compatible columns plus job_id,
                                                                      prep_seed, pose_sha256); failed jobs stay in the table with their status

Identical jobs (same molecule, Vina seed, ligand-preparation seed) are docked once however many sets or scopes contain them, so shared
comparators (reference actives and random molecules) are not recomputed for the second preparation policy.

Pass 1 docks every distinct molecule once (Vina seed 42); pass 2 re-docks the replicate-flagged molecules with two more seeds (engine noise).
Order: molecules of the groups named in --priority first (default M,S,P: the groups with held-out comparators, which H3 needs, come first), each group
in a fixed shuffled order, so an early stop leaves the most important sets complete. The ligand-preparation seed (3-D starting conformer) is --seed for
every job and is part of the job identity.

Resume policy (cdk2moo/dock_jobs.py): finished jobs are skipped; recorded failures are not retried unless --retry-failed; any inconsistency stops the run
unless --repair, which repairs only the REQUESTED jobs (evidence moved to quarantine/, never deleted) and refuses if an unrequested job is inconsistent
(--repair-scope store is the explicit store-wide option). A second concurrent writer is refused.

A Vina score is an estimate of binding energy in kcal/mol, more negative = "better". Its scoring function includes an approximate torsional term but no
rigorous treatment of binding entropy, no explicit water, and the receptor is rigid. It tends to improve with molecule size.

Run (inside a policy scope, after 08c):
    CDK2_CAMPAIGN=<id> CDK2_SCOPE=legacy python scripts/08d_dock_sets.py
    ... python scripts/08d_dock_sets.py --store /tmp/dock_smoke_store --limit 3        # smoke test
"""

import argparse
import functools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cdk2moo import config
from cdk2moo.dock_jobs import ResumeRefused, check_manifest, run_jobs, export, print_progress
from cdk2moo.dock_validate import reconcile, anomalies
from cdk2moo.docking import provenance, prepare_ligand, dock, write_grid_maps, validate_grid_maps


def main():
    config.require_campaign()
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=config.RANDOM_SEED, help="first Vina seed and the ligand-preparation seed")
    ap.add_argument("--store", default=str(config.DOCKING_STORE_DIR or ""), help="docking store directory (default: the campaign's shared store)")
    ap.add_argument("--priority", default="M,S,P", help="group order; molecules of earlier groups are docked first")
    ap.add_argument("--limit", type=int, default=None, help="dock only the first N jobs (smoke test)")
    ap.add_argument("--retry-failed", action="store_true", help="redo requested jobs whose recorded status is a failure")
    ap.add_argument("--repair", action="store_true", help="quarantine and redo requested jobs found inconsistent")
    ap.add_argument("--repair-scope", choices=["requested", "store"], default="requested")
    ap.add_argument("--overwrite-export", action="store_true", help="allow replacing results/08d_export")
    args = ap.parse_args()
    if not args.store:
        sys.exit("no docking store: select a campaign or pass --store")
    store = Path(args.store)
    export_dir = config.RESULTS_DIR / "08d_export"
    if export_dir.exists() and not args.overwrite_export:
        sys.exit(f"{export_dir} exists; this scope's export is not replaced without --overwrite-export")
    seeds = [args.seed, args.seed + 1, args.seed + 2]

    sets = pd.read_csv(config.RESULTS_DIR / "08c_docking_sets.csv")
    group_rank = {g: i for i, g in enumerate(args.priority.split(","))}
    first_group = sets.groupby("smiles")["group"].apply(lambda g: min(group_rank.get(x, 99) for x in g))
    order = []
    for rank in sorted(set(first_group)):
        members = sorted(first_group[first_group == rank].index)
        order += [members[i] for i in np.random.default_rng([args.seed, rank]).permutation(len(members))]
    jobs = [(s, seeds[0]) for s in order]
    replicates = sorted(sets.loc[sets["replicate"], "smiles"].unique())
    jobs += [(s, sd) for s in replicates for sd in seeds[1:]]
    if args.limit:
        jobs = jobs[:args.limit]

    box = config.STRUCTURES_DIR / "receptor_box.txt"
    receptor = config.STRUCTURES_DIR / "receptor.pdbqt"
    from cdk2moo.docking import read_box
    box_values = read_box(box)
    try:
        state = check_manifest(store, provenance(prep_seed=args.seed, receptor=receptor, box=box_values))
        maps_dir = store / "maps"
        if maps_dir.exists():
            validate_grid_maps(maps_dir, receptor, box_values)
        else:
            write_grid_maps(maps_dir, receptor, box_values)
        print(f"store {store}: manifest {state}; grid maps validated; {len(jobs)} dockings requested by this scope "
              f"(exhaustiveness {config.DOCK_EXHAUSTIVENESS}, box {config.DOCK_BOX_SIZE} A)", flush=True)
        start = time.time()
        records, repaired = run_jobs(store, jobs, prepare_ligand, functools.partial(dock, maps_prefix=str(maps_dir / "maps")),
                                     prep_seed=args.seed, retry_failed=args.retry_failed, repair=args.repair,
                                     repair_scope=args.repair_scope, progress=print_progress)
    except ResumeRefused as refusal:
        sys.exit(f"REFUSED: {refusal}")
    if repaired:
        print("repaired (moved to quarantine/ and redone):", {k: len(v) for k, v in repaired.items()})

    report = reconcile(store)
    mine = {s for s, _ in jobs}
    statuses = pd.Series([r["status"] for r in report["records"].values() if r["smiles"] in mine]).value_counts().to_dict()
    print(f"\n{time.time() - start:.0f} s this session. This scope's job records by status: {statuses}")
    print(f"store: complete {len(report['complete'])}, recorded failures {len(report['failed'])}, unresolved inconsistencies: "
          f"{anomalies(report) or 'none'}")
    table = export(store, export_dir, only_smiles=set(sets["smiles"]))
    print(f"exported {len(table)} job rows for this scope to {export_dir}")


if __name__ == "__main__":
    main()

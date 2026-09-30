"""Single source of truth for paths, seeds and constants.

Nothing else in this package hardcodes a path, a seed or a fingerprint
parameter. If you need one, import it from here.

The project root is resolved relative to this file, so the package works from
any working directory and from a cold clone at any location on disk.
"""

from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
# This file lives at <root>/src/cdk2moo/config.py, so the root is three
# parents up.
PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"  # downloaded, never edited by hand
PROCESSED_DIR = DATA_DIR / "processed"  # curated, regenerable
STRUCTURES_DIR = DATA_DIR / "structures"  # PDB / PDBQT files

RESULTS_DIR = PROJECT_ROOT / "results"  # tables, JSON metrics
FIGURES_DIR = PROJECT_ROOT / "figures"  # all output figures
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
NOTEBOOKS_DIR = PROJECT_ROOT / "notebooks"

ALL_OUTPUT_DIRS = (RAW_DIR, PROCESSED_DIR, STRUCTURES_DIR, RESULTS_DIR, FIGURES_DIR)


def ensure_dirs() -> None:
    """Create every output directory if it does not already exist.

    Cheap and idempotent. Call at the top of any script that writes to disk.
    """
    for directory in ALL_OUTPUT_DIRS:
        directory.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Randomness
# --------------------------------------------------------------------------
# Every script exposes --seed and defaults to this value. The seed is logged
# alongside the script's outputs so a figure can always be traced back to the
# run that produced it.
RANDOM_SEED = 42

# Stage 3 runs the GA with multiple seeds so trajectory figures carry error
# bands rather than a single cherry-picked curve.
GA_SEEDS = (42, 43, 44, 45, 46)


# --------------------------------------------------------------------------
# Target and structure
# --------------------------------------------------------------------------
# Human CDK2 (cyclin-dependent kinase 2), ATP-competitive site.
PDB_ID = "4KD1"  # CDK2 + dinaciclib, 1.7 A, monomeric
LIGAND_CODE = "1QK"  # dinaciclib, the co-crystallised ligand in 4KD1

# Fallback structure if redocking 1QK into 4KD1 fails. LS1 is more rigid,
# which usually makes redocking easier.
FALLBACK_PDB_ID = "1KE5"
FALLBACK_LIGAND_CODE = "LS1"

# ChEMBL target for CDK2. Expected to be CHEMBL301, but data.py must verify
# this through the web client rather than trusting the constant blindly.
EXPECTED_CHEMBL_TARGET_ID = "CHEMBL301"

# Dinaciclib is also in ChEMBL. It is held out of surrogate training so it can
# serve as an independent reference point. Its ChEMBL ID is deliberately left
# unset: stage 1 resolves it by structure match and writes it back here rather
# than us hardcoding an ID nobody has checked.
DINACICLIB_CHEMBL_ID = None


# --------------------------------------------------------------------------
# Featurisation
# --------------------------------------------------------------------------
# ECFP4 == Morgan fingerprint with radius 2. 2048 bits is the common default.
MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048

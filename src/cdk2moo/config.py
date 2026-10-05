"""
Central configuration: paths, seeds, and constants.

This is the ONLY place these values are defined. Every other module imports
from here, so that reproducing the project means changing one file.
"""

import os
from pathlib import Path

# --- Paths and campaign selection --------------------------------------
# A pipeline script reads and writes inside ONE selected campaign scope; nothing falls back to the historical
# results. Select with environment variables:
#   CDK2_CAMPAIGN=<id>         a campaign in campaigns/<id>/ (see cdk2moo/campaign.py for its layout)
#   CDK2_CAMPAIGN=historical   explicit reanalysis of the historical artifacts in data/, results/, figures/
#   CDK2_SCOPE=shared|legacy|corrected   which part of the campaign (default shared; legacy and corrected are the two
#                                        candidate-preparation policies and are never pooled)
#   CDK2_PROJECT_ROOT=<dir>    only for tests: a different project root
# With no campaign selected the data paths point at a directory that does not exist, and scripts stop in
# require_campaign() with an explanation. Importing this module never creates directories outside a selected campaign.
PROJECT_ROOT = Path(os.environ.get("CDK2_PROJECT_ROOT") or Path(__file__).resolve().parents[2]).resolve()
CAMPAIGN = os.environ.get("CDK2_CAMPAIGN") or None
SCOPE = os.environ.get("CDK2_SCOPE", "shared")
POLICY_SCOPES = ("legacy", "corrected")
KNOWN_SCOPES = ("shared",) + POLICY_SCOPES

if CAMPAIGN is None:
    _unselected = PROJECT_ROOT / "NO_CAMPAIGN_SELECTED"
    DATA_DIR = RAW_DIR = PROCESSED_DIR = STRUCTURES_DIR = INPUT_STRUCTURES_DIR = _unselected / "data"
    RESULTS_DIR = _unselected / "results"
    FIGURES_DIR = _unselected / "figures"
    CAMPAIGN_ROOT = DOCKING_STORE_DIR = None
    SCOPE_POLICY = None
elif CAMPAIGN == "historical":
    DATA_DIR = PROJECT_ROOT / "data"
    RAW_DIR = DATA_DIR / "raw"
    PROCESSED_DIR = DATA_DIR / "processed"
    STRUCTURES_DIR = INPUT_STRUCTURES_DIR = DATA_DIR / "structures"
    RESULTS_DIR = PROJECT_ROOT / "results"
    FIGURES_DIR = PROJECT_ROOT / "figures"
    CAMPAIGN_ROOT = DOCKING_STORE_DIR = None
    SCOPE_POLICY = None
else:
    if SCOPE not in KNOWN_SCOPES:
        raise ValueError(f"CDK2_SCOPE must be one of {KNOWN_SCOPES}, got {SCOPE!r}")
    CAMPAIGN_ROOT = PROJECT_ROOT / "campaigns" / CAMPAIGN
    _scope_root = CAMPAIGN_ROOT / SCOPE
    DATA_DIR = _scope_root
    RAW_DIR = CAMPAIGN_ROOT / "inputs" / "raw"                 # frozen copy of the cached ChEMBL download
    INPUT_STRUCTURES_DIR = CAMPAIGN_ROOT / "inputs" / "structures"   # frozen 4KD1 structure and crystal ligand
    PROCESSED_DIR = _scope_root / "processed"
    STRUCTURES_DIR = _scope_root / "structures"
    RESULTS_DIR = _scope_root / "results"
    FIGURES_DIR = _scope_root / "figures"
    DOCKING_STORE_DIR = CAMPAIGN_ROOT / "docking_store"        # shared by all scopes: identical jobs are docked once
    SCOPE_POLICY = SCOPE if SCOPE in POLICY_SCOPES else None
    for _d in (PROCESSED_DIR, STRUCTURES_DIR, RESULTS_DIR, FIGURES_DIR):
        _d.mkdir(parents=True, exist_ok=True)


def require_campaign():
    """Stop with an explanation if no campaign is selected; print what was selected otherwise."""
    if CAMPAIGN is None:
        raise SystemExit("No campaign selected. Set CDK2_CAMPAIGN=<id> (and CDK2_SCOPE=shared|legacy|corrected), or "
                         "CDK2_CAMPAIGN=historical to reanalyse the historical artifacts explicitly. There is no default, "
                         "so no script can silently read historical files.")
    print(f"[campaign {CAMPAIGN} | scope {SCOPE if CAMPAIGN != 'historical' else 'n/a'}] results -> {RESULTS_DIR}", flush=True)


def default_policy():
    """The candidate-preparation policy a GA script should use when --policy is not given: the scope's, else legacy."""
    return SCOPE_POLICY or "legacy"


def check_policy(policy):
    """In a policy scope, refuse a different policy (the scope's outputs belong to exactly one policy)."""
    if SCOPE_POLICY is not None and policy != SCOPE_POLICY:
        raise SystemExit(f"scope {SCOPE!r} holds policy {SCOPE_POLICY!r} results only; got --policy {policy!r}")


# --- Reproducibility ---------------------------------------------------
RANDOM_SEED = 42

# --- Target identity ---------------------------------------------------
# Human cyclin-dependent kinase 2.
# UniProt accession is the authoritative identifier: ChEMBL target names are
# inconsistent, accessions are not. Verified against the PDBe entry for CDK2.
TARGET_NAME = "CDK2"
TARGET_UNIPROT = "P24941"
TARGET_ORGANISM = "Homo sapiens"


TARGET_CHEMBL_ID = "CHEMBL301"   # verified via UniProt P24941, 2026-09-30
KEEP_TYPES = ["IC50"]
DROP_CURATOR_FLAGGED = True

# --- Structure ---------------------------------------------------------
# 4KD1: CDK2 + dinaciclib. Verified on RCSB 2026-09-30.
#   X-ray, 1.70 A, R-work 0.193 / R-free 0.232, space group P 21 21 21
#   Single chain A, monomeric (C1), Homo sapiens, UniProt P24941, 0 mutations
#   298 deposited residues, 298 MODELLED -> no missing residues anywhere.
#   Martin et al. (2013) ACS Chem Biol 8:2360, PDB DOI 10.2210/pdb4KD1/pdb
PDB_ID = "4KD1"
LIGAND_CODE = "1QK"      # dinaciclib
LIGAND_FORMULA = "C21 H29 N6 O2"
LIGAND_INCHIKEY = "WBUFFOKXERTKGU-SFHVURJKSA-N"
LIGAND_AUTH_SEQ_ID = 302     # residue number of 1QK in chain A
LIGAND_LABEL_ASYM_ID = "C"   # internal chain id of the ligand instance

# Crystal pose of the ligand as SDF, for the redocking RMSD control.
LIGAND_SDF_URL = (
    "https://models.rcsb.org/v1/4kd1/ligand"
    "?auth_seq_id=302&label_asym_id=C&encoding=sdf&filename=4kd1_C_1QK.sdf"
)

# Non-polymer components to strip during receptor preparation.
#   EDO = 1,2-ethanediol (ethylene glycol), a cryoprotectant, not biology.
STRIP_HETATMS = ["EDO", "HOH"]

PDB_ID_FALLBACK = "1KE5"     # more rigid ligand; use if redocking 4KD1 fails
LIGAND_CODE_FALLBACK = "LS1"

# --- Bioactivity curation ----------------------------------------------
MEASUREMENT_TYPES = ["IC50", "Ki", "Kd", "EC50", "Potency", "Inhibition", "AC50"]

# Decided in stage 1 after looking at the counts. Left explicit rather than
# defaulted, because which types we pool is a judgment call we must defend.
KEEP_TYPES = ["IC50"]
KEEP_UNITS = "nM"
KEEP_RELATION = "="

# Assays run on a mutant CDK2 measure a different protein. Matches either the
# word "mutant", or a residue substitution (e.g. F82H = phenylalanine 82 ->
# histidine) directly after "CDK2". A looser letter-digit-letter pattern was
# tried and false-positived on boilerplate like "Multifuge X1R" and "H2O".
MUTANT_REGEX = r"\bmutant\b|\bCDK2\s+[A-Z]\d{1,3}[A-Z]\b"

# B = binding assay. F (functional) and A (ADMET) are a handful of records.
KEEP_ASSAY_TYPES = ["B"]

# --- Featurization -----------------------------------------------------
FP_RADIUS = 2      # Morgan radius 2 == ECFP4
FP_BITS = 2048

# --- Splits ------------------------------------------------------------
# Fraction of molecules held out as the test set, for both the random and the
# scaffold split. Arbitrary but conventional (80/20); see decisions.md D-12.
TEST_FRAC = 0.2

# Split seeds, fixed BEFORE any surrogate was trained (decisions.md D-16).
# Every split-dependent result is reported as median/IQR over all of them;
# no single seed is ever chosen. The first is the project's default seed.
SPLIT_SEEDS = list(range(42, 52))

# --- Surrogate ---------------------------------------------------------
# Random forest settings, fixed once and NOT tuned (decisions.md D-17): the
# surrogate is an instrument, not a result. max_features=1/3 is Breiman's
# default recommendation for regression forests.
RF_N_TREES = 300
RF_MAX_FEATURES = 0.33

# --- Genetic algorithm (decisions.md D-18 .. D-20) ------------------------
GA_SEEDS = [42, 43, 44, 45, 46]   # 5 runs per arm; fixed in advance
GA_POP_SIZE = 100                 # molecules kept each generation
GA_N_CHILDREN = 100               # new molecules proposed each generation
GA_N_GENERATIONS = 50
GA_MUTATION_PROB = 0.5            # chance a child is mutated after crossover

# PRIMARY size window for generated molecules: the central part of the training
# set's heavy-atom counts, the 5th-95th percentile = 20-39 (decisions.md D-19). The
# training set itself spans 5-78 heavy atoms, so this is NOT its full support: 10%
# of the training molecules lie outside it. The window avoids the sparsely
# populated extremes, where the surrogate has few examples.
GA_MIN_HEAVY_ATOMS = 20
GA_MAX_HEAVY_ATOMS = 39
# RELAXED window, used only for the pilot / stress run (--relaxed). Training
# molecules span 5-78 heavy atoms. Never used for primary Stage 5 conclusions.
GA_RELAXED_MIN_HEAVY_ATOMS = 15
GA_RELAXED_MAX_HEAVY_ATOMS = 50
GA_ALLOWED_ELEMENTS = ["C", "N", "O", "F", "S", "Cl", "Br"]  # that mutation may add

# ONE real surrogate and ONE scrambled surrogate, trained once on all 2,016
# molecules and reused by every GA seed and arm (decisions.md D-18). GA seeds
# then vary only the search, never the instrument.
SURROGATE_SEED = RANDOM_SEED      # random_state of the real forest
SCRAMBLE_SEED = 20260930          # seed of the one-off label shuffle AND of its forest

# Fixed linear rescaling of predicted pActivity to a 0-1 fitness term. The same
# for every arm, so arms differ only in which surrogate they follow.
ACTIVITY_LOW = 4.0
ACTIVITY_HIGH = 9.0

# --- Size-conditioned similarity (decisions.md D-22) ----------------------
# Reference stratum for a generated molecule of h heavy atoms: training
# molecules with h +/- SIZE_STRATUM_HALFWIDTH heavy atoms. If that holds fewer
# than SIZE_STRATUM_MIN_N molecules, the half-width grows by 1 until it does.
SIZE_STRATUM_HALFWIDTH = 2
SIZE_STRATUM_MIN_N = 100

# --- Applicability-domain audit (decisions.md D-25) -----------------------
# Regions are defined on the raw ECFP4 (Morgan radius 2, 2048-bit, binary)
# max-Tanimoto to the training set, using the same edges as the Stage 6 table.
AD_RAW_EDGES = [0.0, 0.4, 0.6, 0.8, 1.01]
AD_RAW_NAMES = ["<0.4", "0.4-0.6", "0.6-0.8", ">=0.8"]
# Regions of the size-conditioned percentile (see similarity.py).
AD_PCT_EDGES = [0.0, 10.0, 50.0, 100.01]
AD_PCT_NAMES = ["<10 (novel for size)", "10-50", ">=50 (typical or better)"]
# A cell of the validation data is only used to estimate error if it holds at
# least this many DISTINCT molecules. Arbitrary, the usual rule of thumb.
AD_MIN_UNIQUE = 30
# "Large error" for the detection test: 1 log unit = a 10-fold misprediction,
# twice the ~0.5 literature inter-lab noise (external figure, not measured here).
AD_LARGE_ERROR = 1.0

# --- Held-out positive controls (decisions.md D-26) -----------------------
# Fixed from the dataset's structure, before any control GA was run.
HOLDOUT_QUANTILE = 0.9           # top decile of measured pActivity (>= 8.05 here)
HOLDOUT_MIN_TOP_PER_FAMILY = 5   # scaffold level: a "potent series" has >= 5 top-decile members
# Generated molecule "approaches" a held-out active at this ECFP4 Tanimoto
# (same boundary as the applicability-domain region A/B split; stage 7).
REDISCOVERY_SIM = 0.6
REDISCOVERY_SIMS = [0.5, 0.6, 0.7, 0.8]   # reported as a sensitivity check

# --- Docking (decisions.md D-29) ---------------------------------------------
DOCK_PH = 7.4                 # protonation of the protein, via OpenMM/pdbfixer defaults
DOCK_BOX_SIZE = 22.0          # Angstrom cube centred on the crystal ligand
DOCK_EXHAUSTIVENESS = 8       # Vina search effort (default); higher = slower, slightly more reproducible
DOCK_N_POSES = 9
REDOCK_MAX_RMSD = 2.0         # Angstrom; above this nothing downstream is trustworthy


def effective_parameters():
    """Every constant that can change a result, for run manifests (so a manifest records the parameters actually in force)."""
    names = ["RANDOM_SEED", "SURROGATE_SEED", "SCRAMBLE_SEED", "SPLIT_SEEDS", "TEST_FRAC", "FP_RADIUS", "FP_BITS", "RF_N_TREES",
             "RF_MAX_FEATURES", "GA_SEEDS", "GA_POP_SIZE", "GA_N_CHILDREN", "GA_N_GENERATIONS", "GA_MUTATION_PROB", "GA_MIN_HEAVY_ATOMS",
             "GA_MAX_HEAVY_ATOMS", "GA_RELAXED_MIN_HEAVY_ATOMS", "GA_RELAXED_MAX_HEAVY_ATOMS", "GA_ALLOWED_ELEMENTS", "ACTIVITY_LOW",
             "ACTIVITY_HIGH", "SIZE_STRATUM_HALFWIDTH", "SIZE_STRATUM_MIN_N", "HOLDOUT_QUANTILE", "HOLDOUT_MIN_TOP_PER_FAMILY",
             "REDISCOVERY_SIM", "REDISCOVERY_SIMS", "AD_RAW_EDGES", "AD_PCT_EDGES", "AD_MIN_UNIQUE", "AD_LARGE_ERROR", "DOCK_PH",
             "DOCK_BOX_SIZE", "DOCK_EXHAUSTIVENESS", "DOCK_N_POSES", "REDOCK_MAX_RMSD", "KEEP_TYPES", "KEEP_UNITS", "KEEP_RELATION",
             "KEEP_ASSAY_TYPES", "MUTANT_REGEX", "DROP_CURATOR_FLAGGED", "TARGET_CHEMBL_ID", "PDB_ID"]
    values = {n: globals()[n] for n in names if n in globals()}
    values.update({"fingerprint": "ECFP4: Morgan radius 2, 2048 bits, binary, chirality NOT included",
                   "activity_scaling": "linear over pActivity [ACTIVITY_LOW, ACTIVITY_HIGH], clipped to 0-1, same for every arm",
                   "standardization": "cleanup, largest fragment, neutralise, canonical tautomer (RDKit TautomerEnumerator defaults: 1000 transforms), organic-element check",
                   "campaign": CAMPAIGN, "scope": SCOPE if CAMPAIGN not in (None, "historical") else None})
    return values

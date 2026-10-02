"""
Central configuration: paths, seeds, and constants.

This is the ONLY place these values are defined. Every other module imports
from here, so that reproducing the project means changing one file.
"""

from pathlib import Path

# --- Paths -------------------------------------------------------------
# Resolved relative to this file, never hardcoded, so the repo works on any
# machine after a clone.
PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"              # downloaded; never edited
PROCESSED_DIR = DATA_DIR / "processed"  # curated; regenerable
STRUCTURES_DIR = DATA_DIR / "structures"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures"

for _d in (RAW_DIR, PROCESSED_DIR, STRUCTURES_DIR, RESULTS_DIR, FIGURES_DIR):
    _d.mkdir(parents=True, exist_ok=True)

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

# PRIMARY size window for generated molecules: the training set's observed
# support (5th-95th percentile of heavy-atom count = 20-39; decisions.md D-19).
# Outside it the surrogate has no training examples of that size, so any
# "improvement" there would be extrapolation by construction.
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

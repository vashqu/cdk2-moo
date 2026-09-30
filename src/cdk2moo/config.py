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

# CDK2 multi-objective molecular optimization

A 4-day proof-of-concept for a PhD interview in computational drug discovery.
The author is a physicist / ML researcher with **no cheminformatics background**.
The project's value is methodological rigour, not chemistry expertise.

---

## The scientific question

> When a molecular optimizer is driven by a machine-learned activity surrogate
> under multi-objective drug-likeness constraints, does improvement in the
> objective reflect progress toward the target, or exploitation of the
> surrogate's extrapolation error?

### Hypotheses

- **H1** Predicted pActivity rises monotonically across GA generations. *(Expected, trivial — the setup, not the finding.)*
- **H2** That rise coincides with increasing distance from the surrogate's training distribution: max Tanimoto to training set *falls* as predicted activity climbs.
- **H3** The gain does not transfer to an orthogonal signal: Vina scores of optimized molecules do not exceed those of held-out known actives.
- **H4** A GA driven by a **label-scrambled** surrogate shows a similar-looking optimization curve — so the curve alone is not evidence of anything.
- **H5** Constraining the GA to stay near the training distribution attenuates H1 but improves transfer in H3. Measure the exchange rate.

**H4 and H5 are the point.** H1 is the setup.

---

## CRITICAL: what counts as success

A negative or null result is a **successful outcome**. This project is designed
to expose a failure mode. Do not "fix" results that look bad.

Never do any of the following without being explicitly asked:

- Change a random seed because a figure looks nicer.
- Drop outliers, filter molecules, or narrow an axis range to clean up a trend.
- Swap a scaffold split for a random split because the metrics are better.
- Tune the surrogate until R² looks impressive. The surrogate is an *instrument*, not a result. A mediocre surrogate is expected and fine.
- Report a metric without its control.
- Re-run with different parameters and present only the final run.

If a result looks wrong, **say so and propose a diagnostic**. Do not silently
make it look right.

Conversely: **actively hunt for leakage and artifacts in our own pipeline.**
Duplicate structures across splits, scaffold overlap, aggregation before
splitting, size confounds, assay heterogeneity. Surface them when noticed.
Three real confounds have already been found this way (see Findings below) —
that hunting is the most valuable thing happening in this project.

---

## Code style — read this before writing any code

**A professor will read this code. Optimise for that.**

- Simple, boring, linear code. No clever abstractions.
- No class hierarchies, no factories, no dependency injection, no plugin systems, no decorators beyond the obvious.
- No config frameworks, no CLI frameworks beyond `argparse`.
- No premature generalisation. If something is used once, write it once, inline.
- Functions do one thing and have a docstring saying what and why.
- Prefer an explicit loop over a dense comprehension when the loop is clearer.
- Comments explain **why**, not what. Chemistry rationale belongs in comments.
- No progress bars, no logging framework, no colour output. `print()` is fine.
- Type hints optional; use them only where they clarify.
- If a file exceeds ~200 lines, it is probably doing too much.

The test: could a computational chemist who does not write Python daily read a
module top to bottom and understand what it does and why? If not, simplify.

---

## Working style

The author is learning this domain. Therefore:

- **Explain chemistry concepts when they first appear.** Assume zero prior cheminformatics. Do not assume familiarity with tautomer, scaffold, PAINS, protonation state, PDBQT, rotatable bond, Murcko scaffold.
- **One stage at a time.** Each stage runs and is inspected before the next is written. Do not generate the whole pipeline in one pass.
- **Flag which choices are standard practice, which are judgment calls, and which are arbitrary.** The author must defend these in an interview.
- When a result is surprising, investigate before proceeding.

---

## decisions.md — maintain this continuously

`decisions.md` at the repo root is a running log of every judgment call. It is
the interview preparation document and it matters as much as the code.

Append an entry whenever a choice is made that could have gone another way.
Format:

```markdown
## D-07: Dropped mutant-CDK2 assay records

**Date:** 2026-09-30
**Stage:** 2 (curation)

**Decision.** Exclude records whose assay_description indicates a mutant CDK2
(F82H, L83V/H84D, analogue-sensitive variants).

**Why.** Hinge residues F82/L83/H84 form the hydrogen bonds an ATP-competitive
inhibitor makes. Mutating them changes the binding site, so these are
measurements against a different protein. One compound spanned 3.8 log units
across mutant and wild-type assays in a single paper.

**Alternatives considered.** Keep them and add a mutant flag as a feature —
rejected, too few records to learn from and it muddles the label.

**Cost.** N records, M compounds.

**How I'd defend it.** [one or two sentences the author can say out loud]

**What would change my mind.** [what evidence would reverse this]
```

Entries are short. The "how I'd defend it" line is the point.

---

## Findings so far (do not rediscover these)

1. **Verified target:** `CHEMBL301`, matched by UniProt **P24941**, `SINGLE PROTEIN`, Homo sapiens. Name matching was rejected as fragile.

2. **pActivity is correct.** Our `9 - log10(nM)` agrees with ChEMBL's `pchembl_value` on all 2,721 comparable records (max diff 0.0055, their rounding). Unit handling is not a risk.

3. **19 ChEMBL IDs collapse to other structures** after standardization (salts, tautomers, charge states). Deduplicating by ID rather than InChIKey would have leaked these across a train/test split.

4. **The noise floor cannot be measured from this data.** 439 compounds have >1 record (443 before the mutant and assay-type filters) but median within-compound SD is exactly 0.000 — most "replicates" are the same measurement curated twice. Even cross-document replicates (268 compounds) show median SD 0.000. Do not quote an empirical noise floor. Cite the literature figure (~0.5 log units inter-lab IC50 reproducibility) as *external*, clearly labelled as not measured here.

5. **Three assay confounds found in `assay_description`:**
   - **Mutant CDK2** (F82H, L83V/H84D, K89T…) — different protein, hinge residues mutated. Up to 3.8 log units within one paper. **Filtered out** (10 records, 4 compounds; all 4 also have wild-type records; D-10).
   - **ATP concentration** (12.5 µM → 8.27, 100 µM → 7.10 for the same compound). Cheng–Prusoff in action. Mostly unstated in descriptions, so it cannot be corrected — **document as an irreducible noise source**, with this example.
   - **Cyclin A present.** Many assays use CDK2/CycA3 despite the SINGLE PROTEIN target filter. CDK2 alone is catalytically near-inactive, so enzymatic assays essentially require the complex. We dock into monomeric CDK2 (4KD1). Cyclin binds away from the hinge so this is tolerable for ATP-competitive binders — **state as a limitation**, do not filter.

6. **`assay_type`**: IC50 set keeps B only (D-11). This was *not* effect-free: F/A removed 35 records and 11 compounds, mixing cell-cytotoxicity assays (rightly dropped) with ~12 real enzyme assays mislabelled F/A. The Ki/Kd external set has **no** assay-type filter, because 404 of its 718 records are one kinome Kd panel (CHEMBL1201862) that ChEMBL labels F.

---

## Target and structure

- **Protein:** human CDK2, ATP-competitive site. UniProt P24941.
- **ChEMBL target:** CHEMBL301.
- **Structure:** PDB **4KD1** — CDK2 + dinaciclib. Verified on RCSB:
  - X-ray, **1.70 Å**, R-work 0.193 / R-free 0.232, P 21 21 21
  - Single chain A, monomeric (C1), 0 mutations
  - **298 deposited residues, 298 modelled — no missing residues anywhere.** No loop rebuilding needed.
  - Ligand `1QK` = dinaciclib, C21 H29 N6 O2, InChIKey WBUFFOKXERTKGU-SFHVURJKSA-N, at chain A residue 302.
  - Deposited as a **1-hydroxypyridinium** — the N-oxide in protonated, formally cationic form. Protonation state is a real decision at docking time; "the crystallographers modelled it as the cation" is the defensible answer.
  - Other heteroatom: `EDO` (ethylene glycol, cryoprotectant). Strip it and waters.
  - Crystal pose available directly as SDF — see `LIGAND_SDF_URL` in config. This is the redocking RMSD reference.
- **Fallback structure** if redocking fails: `1KE5` (ligand `LS1`, more rigid).

---

## Environment

conda env `cdk2`, Python 3.11, macOS **arm64** (MacBook M3). **CPU only.**

conda-forge: `rdkit`, `vina`, `pdbfixer`, `numpy`, `pandas`, `scikit-learn`, `matplotlib`, `seaborn`, `jupyterlab`
pip: `meeko`, `posebusters`, `chembl_webresource_client`

### Hard constraints

- No GPU. No deep generative models. No large model training.
- No molecular dynamics, no free-energy calculations.
- Total compute budget under ~2 hours across the whole project.
- **Never use ADFR Suite or MGLTools** — no arm64 builds exist. Receptor prep goes through `mk_prepare_receptor.py --read_pdb` (Meeko).
- conda-forge packages install before pip packages. **Never `pip install vina`** — it builds from source and needs Boost/SWIG.
- The ChEMBL API drops connections; the raw cache exists so a failure costs one command, not the data.

---

## Repo layout and conventions

```
src/cdk2moo/     importable logic. Module name = the concept.
scripts/         pipeline stages that READ files and WRITE files.
                 NN_verb_noun.py, numbered in execution order.
notebooks/       exploration and inspection. NOT numbered. Nothing load-bearing.
data/raw/        downloaded. Never edited. Cached API responses.
data/processed/  curated. Regenerable by re-running scripts.
data/structures/ PDB, PDBQT, SDF
results/         metrics as JSON/CSV
figures/         every plot
decisions.md     the judgment-call log
```

- A notebook that writes files or feeds a later stage is a script in the wrong place — promote its logic to `src/`.
- A `def` never belongs in `scripts/`.
- Stages communicate **only through files**, never in memory.
- Every script takes `--seed`, defaulting to `config.RANDOM_SEED`. Log it.
- No hardcoded absolute paths. `config.py` is the only place paths, seeds and constants live.
- Config values representing undecided choices default to `None` so an unset decision crashes rather than silently defaulting.

---

## Pipeline stages

| # | Stage | Status |
|---|-------|--------|
| 0 | Environment, structure verification | **done** |
| 1 | ChEMBL raw download (cached) | **done** — 5,367 records, 4,262 compounds |
| 2 | Curation: filter, standardize, aggregate | **done** — 2,016 molecules + 602 Ki/Kd external (mutant filter and assay-type rules applied). |
| 3 | Featurization (ECFP4) + random, scaffold and paper splits, 10 predetermined seeds | **done** (D-13..D-16) |
| 4 | Surrogate: RandomForest; three splits × 10 seeds; scrambled-label, mean and size-only controls | **done** (D-17); median/IQR reported |
| 5 | GB-GA optimization, 4 arms × 5 seeds | **done** (D-18..D-23): primary window 20–39, pilot 15–50 kept; H2 not supported, see D-23; held-out positive controls (D-26) and scaffold-trained robustness across 5 splits (D-27, D-28): 3 of 6 claims stable, 3 fragile |
| 6 | Pareto analysis over (pActivity, QED, SA) | **done** (D-24); overlays are stand-ins, see D-24 |
| 7 | Applicability-domain audit + structural-alert audit | **done** (D-25) |
| 8 | Docking (Vina) — **timeboxed, droppable** | **done** (D-29, D-30; PoseBusters re-audited with all 22 checks in D-34): redocking 0.65 A; Vina does not discriminate actives from decoys, H3 inconclusive |
| 9 | Pose/score audit (PoseBusters, size confound, decoy ROC) | **done** (D-30): PoseBusters 99% pass, size confound quantified, decoy AUC 0.53 |
| 10 | Mitigation arm (H5) | **done** (D-31, D-32): floor is nearly free in predicted activity; H5 not supported as worded |
| 11 | Figures, README, slide deck | figures 1-5 and Table 1 **done** (`figures/final/`); README current; slide deck left to the author |

Stages 1–7 + 10 are a complete, presentable project on their own. Stages 8–9
are upside. **If docking is not working by the Day 3 midpoint, drop it and say
so in the deck.**

---

## Required controls

All of these must exist before the project is done:

1. Held-out known actives (top-decile ChEMBL, excluded from training)
2. Random ChEMBL sample (floor)
3. Property-matched decoys
4. Scaffold split reported alongside random split
5. Scrambled-label surrogate + its GA arm
6. Single-objective ablation arms
7. **Redock the native ligand and report RMSD to the crystal pose** — FIRST step of stage 8. RMSD > 2 Å means nothing downstream is trustworthy and must be reported as such.
8. Heavy-atom-count-only "scoring function" baseline
9. Vina replicate runs with different seeds (engine noise)
10. Multiple GA seeds — error bands on every trajectory figure

---

## Figures (the deliverable)

1. **Surrogate honesty** — predicted vs measured, random split vs scaffold split side by side, plus the scrambled-label control.
2. **Optimization trajectories** — 4 arms, mean ± SD over 5 seeds.
3. **Pareto front** over (pActivity, QED, SA), with known actives and random ChEMBL overlaid. Include rendered structures of 3–4 Pareto-optimal molecules.
4. **Reward hacking** — predicted pActivity vs max Tanimoto to training set, coloured by generation; inset = empirical surrogate error vs similarity on held-out real molecules.
5. **Orthogonal validation** — Vina score distributions per molecule set; Vina vs heavy-atom count; PoseBusters pass rate.

Plus **Table 1**: per-set summary (n, predicted pActivity, Vina, QED, SA, MW,
cLogP, internal diversity, %PAINS, %PoseBusters-pass, max Tanimoto to training).

---

## Limitations to state explicitly

- Assay heterogeneity: ATP concentration varies and is mostly unrecorded; demonstrated at 1.17 log units for one compound in one paper.
- Ligand data is largely against CDK2/cyclin A; docking is into monomeric CDK2.
- Censored `>` records dropped (31% of raw) — the training set is biased toward actives, so the surrogate has little experience of weak compounds.
- No empirical noise floor obtainable from this data.
- Docking scores are not binding free energies: rigid receptor, no explicit water, no rigorous treatment of binding entropy (Vina has only an approximate rotatable-bond penalty), no tautomer/protonation enumeration.
- Single receptor conformation; induced fit ignored.
- SA score is a fragment-frequency heuristic, not a synthetic route.
- QED optimization is partly circular — it rewards resemblance to existing drugs.
- ECFP-based models struggle with activity cliffs: molecules with identical fingerprints (e.g. stereoisomers; 166 curated molecules have an identical-fingerprint twin) cannot be told apart, and similarity-based learners tend to smooth over sharp activity changes between near neighbours. A cliff between molecules with different fingerprints is representable in principle but hard to learn from this much data.
- GB-GA mutations are chemically naive graph edits.
- No ADMET, selectivity or toxicity. No experimental validation of anything.

---

## Non-goals (future work, not in scope)

NSGA-II, uncertainty-aware acquisition, tautomer/protonation enumeration,
ensemble docking, chemical language models with RL, retrosynthesis
(AiZynthFinder), time-split validation, MM-GBSA, FEP, MD.
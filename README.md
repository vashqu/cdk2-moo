# CDK2 multi-objective molecular optimization

A proof-of-concept on human CDK2 (cyclin-dependent kinase 2), ATP-competitive site.

## The scientific question

> When a molecular optimizer is driven by a machine-learned activity surrogate under
> multi-objective drug-likeness constraints, does improvement in the objective reflect
> progress toward the target, or exploitation of the surrogate's extrapolation error?

### Hypotheses

- **H1** Predicted pActivity rises monotonically across GA generations. *(Expected, trivial.)*
- **H2** That rise coincides with increasing distance from the surrogate's training
  distribution — max Tanimoto to the training set *falls* as predicted activity climbs.
- **H3** The gain does not transfer to an orthogonal signal: Vina scores of optimized
  molecules do not exceed those of held-out known actives.
- **H4** A GA driven by a **label-scrambled** surrogate shows a similar-looking
  optimization curve — so the curve alone is not evidence of anything.
- **H5** Constraining the GA to stay near the training distribution attenuates H1 but
  improves transfer in H3. Measure the exchange rate.

H4 and H5 are the point. H1 is the setup. **A null or negative result is a successful
outcome here** — the project is built to expose a failure mode, not to beat a benchmark.

## Pipeline stages

| # | Stage | One line |
|---|-------|----------|
| 0 | Environment + structure inspection | Check the env, pull PDB `4KD1`, look at the binding site and ligand `1QK`. |
| 1 | ChEMBL download, curation, standardization | Fetch CDK2 bioactivities, strip salts, neutralize, canonicalize tautomers, de-duplicate, aggregate replicates, convert to pActivity. |
| 2 | Surrogate | ECFP4 → RandomForest, evaluated under random *and* scaffold splits, plus a label-scrambled control. |
| 3 | GB-GA optimization | Graph-based genetic algorithm, 4 arms × 5 seeds, logging every molecule per generation. |
| 4 | Pareto analysis | Non-dominated fronts over (pActivity, QED, SA), with known actives and random ChEMBL overlaid. |
| 5 | Applicability-domain audit | Max Tanimoto to training vs predicted activity; empirical surrogate error vs similarity on held-out real molecules. |
| 6 | Docking (Vina) | **Timeboxed and droppable.** Native-ligand redocking first; if RMSD > 2 Å, downstream docking numbers are reported as untrusted. |
| 7 | Pose/score audit | PoseBusters pass rate, heavy-atom-count confound, ROC vs property-matched decoys. |
| 8 | Mitigation arm (H5) | Re-run the GA under a similarity constraint and measure what the constraint costs and buys. |

Stages 1–5 and 8 are a complete project on their own. Stages 6–7 are upside.

## Reproducing

```bash
conda env create -f environment.yml
conda activate cdk2
pip install -e .
```

Then run everything in `scripts/` in numeric order:

```bash
for s in scripts/[0-9][0-9]_*.py; do python "$s"; done
```

Every script takes `--seed` (defaulting to `cdk2moo.config.RANDOM_SEED`), logs the seed
it used, and writes its outputs to disk. Nothing is passed between stages in memory, so
any stage can be re-run on its own.

## Layout

```
src/cdk2moo/   importable modules; config.py is the only place paths,
               seeds and constants are defined
scripts/       00_, 01_, ... runnable in numeric order
data/raw/      downloaded, never edited (git-ignored)
data/processed/curated, regenerable (git-ignored)
data/structures/PDB files
results/       tables, JSON metrics (git-ignored)
figures/       output figures (git-ignored)
notebooks/     exploration only, nothing load-bearing
```

`src/cdk2moo/docking.py` is deliberately isolated: nothing outside it and the docking
scripts imports from it, so stage 6 can be deleted wholesale without breaking the rest.

## Constraints

CPU only (macOS arm64), no GPU, no deep generative models, no MD or free-energy
calculations, total compute budget under ~2 hours. Receptor prep goes through Meeko's
`mk_prepare_receptor.py --read_pdb`; ADFR Suite and MGLTools are never used (no arm64
builds).

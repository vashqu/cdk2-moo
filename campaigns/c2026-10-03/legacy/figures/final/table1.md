| set | n | n with Vina | pActivity | basis | Vina | QED | SA | MW | cLogP | diversity | % PAINS | % PB pass | max Tanimoto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| P:training_top_decile | 100 | 100 | 8.40 | measured | -8.31 | 0.53 | 2.78 | 405 | 3.35 | 0.83 | 15 | 99 | 0.78 |
| M:heldout_actives | 100 | 100 | 8.40 | measured | -8.17 | 0.52 | 2.91 | 410 | 3.44 | 0.83 | 9 | 99 | 0.73 |
| S:heldout_actives | 93 | 93 | 8.40 | measured | -8.18 | 0.54 | 3.03 | 428 | 3.38 | 0.75 | 9 | 100 | 0.67 |
| P:random_training | 100 | 100 | 6.59 | measured | -8.23 | 0.52 | 2.86 | 428 | 3.87 | 0.86 | 13 | 99 | 0.78 |
| M:random_remaining | 100 | 100 | 6.28 | measured | -8.28 | 0.55 | 2.93 | 417 | 3.45 | 0.87 | 7 | 96 | 0.78 |
| S:random_remaining | 100 | 100 | 6.01 | measured | -8.44 | 0.54 | 2.78 | 418 | 3.63 | 0.87 | 7 | 99 | 0.79 |
| P:multi_real | 100 | 100 | 7.70 | predicted | -7.78 | 0.90 | 2.06 | 305 | 2.74 | 0.72 | 1 | 99 | 0.38 |
| P:multi_scrambled | 100 | 100 | 6.16 | predicted | -7.64 | 0.91 | 2.22 | 317 | 3.03 | 0.79 | 6 | 99 | 0.37 |
| P:activity_only | 100 | 100 | 8.58 | predicted | -8.86 | 0.30 | 3.06 | 484 | 5.01 | 0.50 | 10 | 97 | 0.71 |
| P:druglike_only | 100 | 100 | 6.06 | predicted | -7.54 | 0.93 | 1.72 | 298 | 3.15 | 0.78 | 11 | 99 | 0.38 |
| P:pareto_front | 124 | 124 | 7.59 | predicted | -7.88 | 0.92 | 1.85 | 306 | 3.05 | 0.78 | 6 | 99 | 0.43 |
| P:decoys | 100 | 99 | 6.37 | predicted | -8.11 | 0.62 | 2.95 | 396 | 3.24 | 0.86 | 7 | 94 | 0.48 |
| M:multi_real | 100 | 100 | 7.22 | predicted | -7.97 | 0.91 | 1.89 | 304 | 2.68 | 0.67 | 2 | 97 | 0.35 |
| M:multi_scrambled | 100 | 100 | 6.05 | predicted | -7.80 | 0.92 | 1.80 | 311 | 2.74 | 0.74 | 7 | 99 | 0.37 |
| M:activity_only | 100 | 100 | 7.84 | predicted | -8.16 | 0.29 | 4.09 | 522 | 3.56 | 0.60 | 9 | 89 | 0.70 |
| M:druglike_only | 100 | 100 | 6.06 | predicted | -7.55 | 0.93 | 1.72 | 298 | 2.94 | 0.76 | 23 | 100 | 0.36 |
| M:decoys | 100 | 100 | 6.25 | predicted | -8.02 | 0.63 | 2.86 | 393 | 3.28 | 0.86 | 6 | 100 | 0.48 |
| S:multi_real | 100 | 100 | 7.43 | predicted | -7.94 | 0.91 | 2.35 | 319 | 2.66 | 0.79 | 6 | 99 | 0.42 |
| S:multi_scrambled | 100 | 100 | 5.93 | predicted | -7.38 | 0.92 | 2.12 | 308 | 2.92 | 0.82 | 15 | 99 | 0.32 |
| S:activity_only | 100 | 100 | 8.39 | predicted | -8.53 | 0.19 | 3.31 | 506 | 3.97 | 0.68 | 84 | 98 | 0.68 |
| S:druglike_only | 100 | 100 | 6.09 | predicted | -7.70 | 0.93 | 1.77 | 297 | 3.18 | 0.78 | 17 | 99 | 0.35 |
| S:decoys | 93 | 93 | 6.42 | predicted | -8.08 | 0.62 | 2.90 | 400 | 3.34 | 0.86 | 5 | 98 | 0.49 |

Table 1. Campaign c2026-10-03, legacy policy scope. Per-set summary of docked molecules (medians; diversity = mean pairwise 1 - ECFP4 Tanimoto). Groups: P = surrogate trained on all 2,016 molecules (its comparator actives were in training); M = molecule-level holdout; S = scaffold-level holdout. pActivity is measured for sets drawn from the curated data and predicted by the group's surrogate otherwise; not comparable. Vina in kcal/mol (more negative = better), median over the n with Vina molecules. % PB pass = passes all 22 PoseBusters checks / n in set. Max Tanimoto is raw ECFP4 (Morgan r=2, 2048-bit) to the group's training molecules, excluding the molecule itself for sets drawn from training. random_training / random_remaining are stand-ins for a random ChEMBL sample, not a ChEMBL sample.

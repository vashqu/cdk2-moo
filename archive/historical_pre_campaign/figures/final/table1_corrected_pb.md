| set | n | pActivity | basis | Vina | QED | SA | MW | cLogP | diversity | % PAINS | % PB pass | max Tanimoto |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| known_actives | 100 | 8.40 | measured | -8.15 | 0.52 | 2.94 | 417 | 3.25 | 0.82 | 9 | 100 | 0.78 |
| random_training | 100 | 6.51 | measured | -8.43 | 0.51 | 2.83 | 431 | 3.70 | 0.87 | 5 | 100 | 0.80 |
| decoys | 100 | 6.44 | predicted | -8.08 | 0.61 | 2.89 | 395 | 3.10 | 0.87 | 10 | 98 | 0.50 |
| multi_real | 100 | 7.74 | predicted | -7.82 | 0.91 | 2.26 | 304 | 2.84 | 0.73 | 1 | 98 | 0.40 |
| multi_scrambled | 100 | 5.93 | predicted | -7.70 | 0.93 | 2.14 | 314 | 2.93 | 0.81 | 9 | 100 | 0.35 |
| activity_only | 100 | 8.59 | predicted | -8.82 | 0.32 | 3.14 | 479 | 5.08 | 0.51 | 6 | 92 | 0.70 |
| druglike_only | 100 | 6.00 | predicted | -7.64 | 0.93 | 1.73 | 295 | 3.19 | 0.77 | 22 | 100 | 0.36 |
| pareto_front | 101 | 7.85 | predicted | -8.08 | 0.91 | 2.23 | 308 | 3.07 | 0.76 | 5 | 98 | 0.43 |

Table 1. Per-set summary of docked molecules (medians; diversity = mean pairwise 1 - ECFP4 Tanimoto). pActivity is measured for known_actives and random_training and predicted by the all-data random forest for the other sets; not comparable. Vina in kcal/mol (more negative = better). Max Tanimoto is the raw ECFP4 (Morgan r=2, 2048-bit) similarity to the 2,016 training molecules, leave-one-out for training members. random_training is a stand-in for a random ChEMBL sample, not a ChEMBL sample.

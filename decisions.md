# Decisions log

Every judgment call in the project, with the reasoning and a line to say out
loud in an interview. Entries are appended as choices are made.

Each entry is tagged **[standard]**, **[judgment]** or **[arbitrary]**:
standard = what most people in the field would do; judgment = defensible either
way, must be argued; arbitrary = had to pick something, no real argument.

---

## D-01: Target identified by UniProt accession, not by name

**Date:** 2026-09-30
**Stage:** 1 (download) — **[standard]**

**Decision.** The ChEMBL target is `CHEMBL301`, chosen because its protein
component has UniProt accession P24941, its type is `SINGLE PROTEIN`, and its
organism is *Homo sapiens*.

**Why.** A ChEMBL name search for "CDK2" returns protein complexes
(CDK2/cyclin A), non-human orthologues and cell lines. Names are free text and
inconsistent; a UniProt accession is a stable identifier for one specific
protein sequence.

**Alternatives considered.** Take the top hit of the name search — rejected as
fragile; the ranking is not a statement about identity.

**Cost.** None in data terms. One-off effort in `notebooks/explore_target.ipynb`.

**How I'd defend it.** "I matched on the protein's sequence identifier, not its
label. A name search returns complexes and other species, and I didn't want my
dataset's identity to depend on a search ranking."

**What would change my mind.** Nothing realistic; if ChEMBL re-mapped CHEMBL301
to a different accession, the check in the notebook would show it.

---

## D-02: SINGLE PROTEIN targets only (CDK2/cyclin A complex excluded as a target)

**Date:** 2026-09-30
**Stage:** 1 (download) — **[judgment]**

**Decision.** Download only records filed under the `SINGLE PROTEIN` target
CHEMBL301, not under the separate CDK2/cyclin A *protein complex* target.

**Why.** It keeps the target definition unambiguous and identical to the
protein we dock into. Records filed against the complex are, by ChEMBL's own
labelling, measurements of a different biological object.

**Alternatives considered.** Pool the complex-target records too — more data,
but the target identity becomes a mixture, and it buys little because the
single-protein records already include many CDK2/cyclin A assays (see below).

**Cost.** Not measured: we never downloaded the complex target's records, so we
do not know how many compounds we forgo. Could be counted with one API call.

**Important caveat (a limitation, not a filter).** The SINGLE PROTEIN filter
does **not** remove cyclin A assays. Many records filed under CHEMBL301 were
actually run on CDK2/cyclin A, because CDK2 alone is catalytically near-inactive.
So "single protein" describes ChEMBL's annotation, not the assay. We dock into
monomeric CDK2 (4KD1); cyclin binds away from the hinge, so this is tolerable
for ATP-competitive binders. Stated as a limitation.

**How I'd defend it.** "I used the single-protein target so the label matches
the structure I dock into. I know many of those assays still contain cyclin A;
I've listed that as a limitation rather than pretend the filter removed it."

**What would change my mind.** If the dataset were too small to train on, I
would add the complex records with a flag and test whether they shift labels.

---

## D-03: IC50 only — not pooled with Ki/Kd

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[standard]**

**Decision.** The modelling set uses `standard_type == IC50` only. Ki and Kd
records are curated separately and used only as an external check, restricted
to compounds with no IC50.

**Why.** IC50 is the concentration that halves enzyme activity *in that
assay*, so it depends on the assay conditions, notably ATP concentration
(Cheng–Prusoff: Ki = IC50 / (1 + [ATP]/Km)). Ki is the assay-independent
binding constant. For an ATP-competitive inhibitor the same compound has a
numerically *larger* IC50 than Ki, by a factor depending on [ATP]. Pooling
would put two differently-offset quantities on one axis.

**Alternatives considered.** Pool everything with a type flag; convert IC50 to
Ki with Cheng–Prusoff — rejected, [ATP] is mostly unrecorded.

**Cost.** Raw funnel: 5,367 records → 3,692 IC50 (1,675 Ki/Kd records set aside).
The Ki/Kd-only compounds are kept as an external set (595 before the mutant
filter).

**How I'd defend it.** "IC50 and Ki aren't the same number for the same
compound — the offset depends on ATP concentration. I kept the label
homogeneous in type and used Ki/Kd as an independent check instead."

**What would change my mind.** If the surrogate were badly data-starved, I'd
test pooling with a type indicator feature and see whether the scaffold-split
error improved.

---

## D-04: Censored ('>') records dropped

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[standard, with a known cost]**

**Decision.** Keep `standard_relation == '='` only. Records like
"IC50 > 10000 nM" are dropped.

**Why.** A regression label needs a number. A censored record says "weaker than
X" — real information but not a value.

**Alternatives considered.** Impute at the threshold; use a censored-regression
loss; classification on a threshold — rejected as out of scope for a
proof-of-concept.

**Cost.** Of the 3,692 IC50 records, 866 are not `=` (23%). (CLAUDE.md's "31% of
raw" is across IC50+Ki+Kd.) The surviving data is biased toward active
compounds: **the surrogate has seen few weak molecules**, which matters for
H2 — it will extrapolate badly toward the inactive side, and the GA may find
"high predicted activity" regions the model has no evidence about.

**How I'd defend it.** "I dropped them because I had no number to regress on,
and I'm explicit that this makes the training set optimistic. It is a named
limitation, not a hidden one."

**What would change my mind.** If the applicability-domain audit (stage 7)
shows the worst errors sit exactly where censored data would have informed the
model.

---

## D-05: Curator-flagged records dropped

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[standard]**

**Decision.** Drop any record with a non-empty `data_validity_comment`
(e.g. "Potential transcription error", "Outside typical range").

**Why.** These are ChEMBL curators' own warnings that the value may be wrong.
Training on records the data source itself distrusts is hard to defend.

**Alternatives considered.** Keep and flag — rejected; too few to learn from.

**Cost.** 73 records (2,798 → 2,725 after the units filter).

**How I'd defend it.** "ChEMBL's curators flagged these as suspect; I took
their word rather than train on values the source doesn't trust."

**What would change my mind.** If the flagged records turned out to be mostly
genuine extreme values (e.g. very potent compounds) — dropping them would
narrow the top of the range, which is precisely where the GA operates.
Worth a one-off look, not yet done.

---

## D-06: Deduplicate by InChIKey of the standardized structure, not by ChEMBL ID

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[standard]**

**Decision.** After standardization (salt strip, neutralise, canonical
tautomer), group records by InChIKey. One row per unique structure.

**Why.** ChEMBL assigns different IDs to a free base and its salt, or to
different tautomer drawings. A random split by ID would put one in train and
one in test: test leakage.

**Alternatives considered.** Deduplicate by `molecule_chembl_id` — rejected;
**19 IDs were found to collapse onto other structures** this way.

**Cost.** 2,721 records → 2,027 unique structures (before the mutant filter).

**How I'd defend it.** "Same molecule, different ID is the textbook source of
train/test leakage. I found 19 instances in this data, so it isn't hypothetical."

**What would change my mind.** InChIKey ignoring a chemically meaningful
distinction (e.g. certain tautomer-sensitive stereochemistry) that we
need to keep. Not seen.

---

## D-07: Replicates aggregated as the median of pActivity

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[judgment]**

**Decision.** A compound measured several times gets one label: the median of
its pActivity values (pActivity = 9 − log10 nM).

**Why.** (i) Median over mean: one measurement at an unusual ATP concentration
can sit a log unit from the rest and drag a mean. (ii) pActivity over raw nM:
averaging concentrations linearly lets the weakest measurement dominate; the
log scale is proportional to binding free energy and makes errors
roughly symmetric.

**Alternatives considered.** Mean of pActivity; geometric mean of nM (this is
the mean of pActivity, same thing); keep every record as a separate row —
rejected, would weight well-studied compounds more heavily and put replicates
on both sides of a split.

**Cost.** Information on spread is kept (`pactivity_std`, `pactivity_range`)
but collapsed out of the label.

**How I'd defend it.** "Median for robustness to the odd outlier assay, and in
log space because that's where the error is roughly symmetric and where the
number means something energetically."

**What would change my mind.** A large, consistent difference between median
and mean labels would suggest outliers matter; I'd examine them rather than
choose the nicer one.

---

## D-08: PDB 4KD1 as the docking structure

**Date:** 2026-09-30
**Stage:** 0 (structure) — **[judgment]**

**Decision.** Dock into 4KD1 (CDK2 + dinaciclib), with 1KE5 as fallback.

**Why.** Verified on RCSB: X-ray at 1.70 Å, R-free 0.232, single monomeric
chain, 0 mutations, 298/298 residues modelled so no loop rebuilding is needed.
The ligand is a potent ATP-site inhibitor, so the pocket is in a
ligand-bound conformation. The crystal pose is available as SDF, giving the
redocking control (required control 7).

**Alternatives considered.** 1KE5 (more rigid ligand) — held as a fallback if
redocking 4KD1 gives RMSD > 2 Å.

**Cost.** One rigid receptor conformation; induced fit ignored. It lacks
cyclin A, while much of the ligand data was measured with it. The ligand is
deposited as a 1-hydroxypyridinium (cationic N-oxide), so protonation state is
a decision at docking time.

**How I'd defend it.** "High resolution, complete, monomeric, unmutated, and the
co-crystal pose lets me test my docking setup by redocking before I trust any
score. Single-conformation docking is a limitation I state."

**What would change my mind.** Redocking RMSD > 2 Å → switch to 1KE5.

---

## D-09: No empirical noise floor is claimed

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[judgment]**

**Decision.** We do not quote a measured label-noise floor. We cite the
literature figure (~0.5 log units inter-lab IC50 reproducibility) as
*external*, labelled as not measured here.

**Why.** 443 compounds have more than one record, but the median within-compound
SD is 1.9e-5 — effectively exactly zero. Most "replicates" are the same
measurement curated twice (the same number entered under two records). Even
cross-document replicates (268 compounds) show median SD 0.000. The real
disagreement is hidden, so a floor computed from this data would be falsely
optimistic, and the surrogate would be judged against a bar that is too low.

**Alternatives considered.** Report the SD anyway — rejected, misleading.
Restrict to compounds with >1 *distinct* value — biased toward disagreement,
no honest denominator.

**Cost.** We cannot say how much of the surrogate's error is irreducible noise.

**How I'd defend it.** "I looked for the noise floor and found the replicates
aren't independent, so I didn't report a number I couldn't stand behind. I cite
the literature value and say plainly that it isn't measured here."

**What would change my mind.** A dataset with genuinely independent replicates
(different labs, raw records).

---

## D-10: Mutant-CDK2 assay records dropped

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[standard]**

**Decision.** Drop records whose `assay_description` matches
`\bmutant\b|\bCDK2\s+[A-Z]\d{1,3}[A-Z]\b` (case-insensitive): the word
"mutant", or a residue substitution such as F82H written directly after
"CDK2". All mutants are dropped, including the two cyclin-interface mutants
(R150A, Y180A), which were already removed by the earlier filters.

**Why.** A mutant CDK2 is a different protein. F82/L83/H84 form the hinge, the
backbone segment whose hydrogen bonds anchor an ATP-competitive inhibitor, and
F80 is the gatekeeper residue at the back of the pocket (F80G is the
analog-sensitive variant). Changing them changes the binding site itself. One
compound spans 3.4 log units (was 3.8 before the filter) across mutant and
wild-type assays.

**Alternatives considered.** Keep them with a mutant flag feature — rejected,
too few records. A looser letter-digit-letter pattern — rejected after
inspection: it false-positived on "Multifuge X1R" and "H2O" in two long
boilerplate descriptions (~294 records). Keyword-only searches (mutation,
variant, gatekeeper, knock…) found nothing the tight pattern missed.

**Cost.** 10 records, 4 compounds. All 4 compounds also have wild-type
records, so no molecule was lost; their labels are now medians over fewer
records. Effect on the dataset is tiny.

**How I'd defend it.** "Those measurements are against a protein whose
inhibitor-binding residues were deliberately changed, so they aren't
wild-type CDK2 data. I checked the regex by hand against every match and
against keywords it might miss."

**What would change my mind.** A mutant description the regex misses (it only
catches mutants that say "mutant" or name the substitution after "CDK2").

---

## D-11: Assay-type filter: B only for the IC50 set, none for the Ki/Kd external set

**Date:** 2026-09-30
**Stage:** 2 (curation) — **[judgment]**

**Decision.** The IC50 modelling set keeps `assay_type == 'B'` only. The Ki/Kd
external set has **no** assay-type filter.

**Why.** ChEMBL's B = biochemical binding/enzyme assay, F = functional
(typically cell-based), A = ADMET. For IC50, the F/A records include
cell-cytotoxicity MTT assays that merely mention CDK2; cell death is not CDK2
binding, so a homogeneous label needs them out. For Ki/Kd, 404 of 718 records
are one document (CHEMBL1201862, "Navigating the Kinome", a 395-compound
kinase panel reporting Kd). ChEMBL labels it F, but Kd is a direct binding
readout and, unlike IC50, does not depend on ATP concentration. The label is
a ChEMBL annotation quirk, not a sign the data is unsuitable.

**Alternatives considered.** One rule for both sets (external shrinks
595 → 227); keep B and F for the IC50 set too — rejected, it lets cell assays in.

**Cost.** IC50 set: 35 records, 13 compounds, 11 lost entirely
(2,027 → 2,016 molecules). This was *not* "for cleanliness, not effect" as
CLAUDE.md assumed. About 12 of the dropped records are genuine enzyme assays
(biotinylated-peptide substrate, Hotspot) that ChEMBL labelled F/A; they are
lost. The two sets are filtered on different rules, which must be stated.
External set: 602 molecules. (Up from 595 before D-10/D-11 because 7 of the 11
compounds that left the IC50 set have Ki/Kd records, so they now qualify as
"no IC50".) The kinome panel is mostly promiscuous kinase inhibitors, so the
external set is chemically unlike the CDK2 training set — useful as an
extrapolation test, weak as a like-for-like check.

**How I'd defend it.** "For IC50 I wanted enzyme assays only, since cell
readouts aren't CDK2 binding. For the external set the F label on the kinome
panel is ChEMBL's annotation, the readout is a binding constant, so I kept it
rather than lose 60% of the set."

**What would change my mind.** Reading the dropped IC50 enzyme-assay
descriptions and finding them equivalent to the B assays (then I'd rescue them
by description, not by type); or finding the kinome panel's Kd values
inconsistent with other Kd data for the same compounds.

---

## D-12: ECFP4, 2048 bits, as the molecular representation

**Date:** 2026-09-30
**Stage:** 3 (features) — **[standard; bit count arbitrary]**

**Decision.** Morgan fingerprint, radius 2 (= ECFP4), 2048 bits, no chirality,
binary (on/off) rather than counts.

**Why.** It is the default representation for classical QSAR and fast, which
matters on a CPU-only budget. Radius 2 is the conventional compromise between
too-local (radius 1) and too-specific (radius 3+) substructures. The
GA later mutates molecules, and fingerprints can be computed for any valid
molecule with no training.

**Alternatives considered.** RDKit descriptors; count fingerprints; learned
graph features — out of scope (no GPU, no deep models).

**Cost.** Bit collisions (two substructures hashed to one bit); no stereochemistry,
so stereoisomers are indistinguishable; activity cliffs — pairs of near-identical
molecules with very different activity — cannot be represented. The 2048
choice is arbitrary.

**How I'd defend it.** "It's the standard baseline, cheap, and its limitations
(no stereochemistry, blind to cliffs) are exactly the kind of blind spots an
optimizer can exploit, which is part of what I'm studying."

**What would change my mind.** Nothing within scope; I'd revisit only to test
whether the conclusions depend on the representation.

---

## D-13: Random and balanced-scaffold splits, 80/20, same seed

**Date:** 2026-09-30
**Stage:** 3 (splits) — **[standard; balanced variant is a judgment]**

**Decision.** Two train/test splits, 80/20, seed 42, saved to
`data/processed/cdk2_splits.csv` and reused by every later stage:
(a) random; (b) scaffold split with whole Bemis-Murcko scaffold groups on
one side. The scaffold split is the "balanced" variant: groups larger than half
the test set go to train, the rest are shuffled by the seed.

**Why.** The random split measures interpolation among relatives; the scaffold
split measures prediction for new ring systems, closer to what an optimizer
asks. Reporting both is a required control. The deterministic textbook version
(largest groups to train, singletons to test) would make the test set all
singletons and ignore the seed.

**Alternatives considered.** Cluster split on Tanimoto similarity; split by
source document; time split (non-goal). 80/20 with no separate validation set:
hyperparameter choices, if any, use cross-validation inside train.

**Cost.** Test scaffolds are not literally independent (see the finding below):
the scaffold split removes identical ring skeletons, not near-identical ones.
Test size is approximate (here exactly 403 = 20.0%).

**Result at seed 42 (test vs train).**

| | random | scaffold |
|---|---|---|
| test molecules whose scaffold is in train | 280 / 403 | 0 / 403 |
| median max Tanimoto to train | 0.787 | 0.639 |
| test molecules with max Tanimoto > 0.7 | 76% | 38% |
| test molecules sharing a source paper with train | 383 / 403 | 376 / 403 |

**How I'd defend it.** "I report both because the gap between them is the
honest measure of how much of the surrogate's accuracy is memorising
analogue series. I verified the scaffold split leaks no scaffold, and I checked
similarity and shared-paper overlap rather than trusting the label 'scaffold
split'."

**What would change my mind.** If a document-grouped split shows a much larger
gap than the scaffold split, it becomes the primary hard split.

---

## D-14: Third split grouped by source paper (follow-up to D-13)

**Date:** 2026-09-30
**Stage:** 3 (splits) — **[judgment]**

**Decision.** Add `paper_split`: molecules grouped by `primary_document` (the
first paper that measured them), same balanced group algorithm, 80/20, seed 42.
Stage 2 now also records `all_documents` (every paper per molecule) so leakage
through secondary papers can be measured. Labels are unchanged; the new column
is the only difference in the curated file.

**Why.** D-13 found 376/403 scaffold-split test molecules came from a paper that
also supplied training molecules. Papers publish analogue series, so the
paper split is the strictest test of extrapolation available from this data.

**Alternatives considered.** Connected components over papers linked by shared
compounds — rejected: reference inhibitors are re-measured everywhere, and one
component holds ~1,128 of the ~2,676 records (42%), too lumpy for an 80/20 split.

**Result at seed 42 (test vs train).**

| | random | scaffold | paper |
|---|---|---|---|
| test molecules with scaffold in train | 280 | 0 | 112 |
| test molecules whose first paper is in train | 383 | 376 | 0 |
| ...whose ANY paper is in train (strict) | 386 | 383 | **166** |
| median max Tanimoto to train | 0.787 | 0.639 | 0.432 |
| test with max Tanimoto > 0.7 | 76% | 38% | 11% |
| pActivity mean train / test | 6.58 / 6.60 | 6.61 / 6.47 | **6.50 / 6.89** |

**Cost / caveats.**
- The paper split still leaks: 166/403 test molecules were also measured in some
  paper that supplied training molecules (267 molecules appear in >1 paper).
- **The test labels are shifted:** mean pActivity 6.89 vs 6.50 in train and a
  wider spread. With lumpy groups, the seed decides which papers land in test.
  Any "surrogate does worse on the paper split" conclusion therefore mixes
  novelty with a different label distribution. Not reseeded (would be
  seed-shopping); to be handled by reporting across several seeds if the paper
  split is used as a headline result.
- Scaffold and paper splits both leave 'big group' rule inactive here: the
  largest group (74 molecules / 131 papers) is below half the test size (201).

**How I'd defend it.** "Scaffold splits are the standard, but I measured that
they still leave near-analogues across the divide, so I added a stricter
paper split and report how much leakage is left in each, rather than assume."

**What would change my mind.** Surrogate metrics that barely differ between
scaffold and paper splits (then the extra split adds little), or a split-seed
sweep showing the paper-split label shift is an artifact of one draw.

---

## D-15: Roles of the three splits

**Date:** 2026-09-30
**Stage:** 3–4 — **[judgment]**

**Decision.** All three splits are kept and reported. **Scaffold = primary
controlled test of chemical extrapolation** (new ring skeletons, same
distribution of assays and labels as far as possible). **Paper = stricter
domain-shift stress test** (new series, new lab, new assay conditions).
**Random = optimistic interpolation baseline.** The applicability-domain audit
(stage 7) analyses error against similarity to the training set *across all
three*, using `results/04_test_predictions.csv`.

**Why.** The scaffold split changes one thing (ring skeleton). The paper split
changes chemistry *and* assay conditions (notably ATP concentration, which is
mostly unrecorded; see CLAUDE.md finding 5), so a drop on it is not purely a
chemistry effect, and its label shift varies widely with the draw
(pActivity test−train from −0.43 to +0.39 across seeds).

**Cost.** Two hard splits to explain. The paper split still leaks (median 113
of 403 test molecules share *some* paper with train, via reference compounds).

**How I'd defend it.** "Each split changes a different thing, so the three
together separate interpolation, chemical novelty, and domain shift. I call
the paper split a stress test, not a cleaner scaffold split, because it also
moves the assay conditions."

**What would change my mind.** Scaffold and paper results being
indistinguishable across seeds.

---

## D-16: Ten predetermined split seeds; report median/IQR, never pick a seed

**Date:** 2026-09-30
**Stage:** 3–4 — **[arbitrary number, principled protocol]**

**Decision.** `config.SPLIT_SEEDS = 42..51`, fixed before any surrogate was
trained (only split statistics had been seen, no model performance). Every
split-dependent metric is reported as median [25th, 75th percentile] over all
ten, for all three splits. `cdk2_splits.csv` holds all ten, in long format.

**Why.** With lumpy group splits, one seed's test set can be unrepresentative:
seed 42's paper split happens to be among the most label-shifted draws
(+0.39 vs median −0.06). Choosing a seed is seed-shopping; a pre-declared list
removes the choice. Ten is enough for a stable IQR and cheap.

**Cost.** Ten times the compute of one split (stage 4: ~4 min, 10 seeds × 3
splits × 4 models). Test sets overlap across seeds, so the ten results are not
independent draws: the IQR understates true uncertainty about "new data".

**How I'd defend it.** "I fixed the seeds in advance and report the spread,
because a single group split can land on an easy or hard draw."

**What would change my mind.** A need for formal confidence intervals would
call for nested/repeated CV with proper independence.

---

## D-17: Surrogate = random forest on ECFP4, fixed hyperparameters, four-way control

**Date:** 2026-09-30
**Stage:** 4 — **[standard model, arbitrary settings]**

**Decision.** `RandomForestRegressor`, 300 trees, `max_features=0.33`, rest at
scikit-learn defaults; **no tuning**. Every test set is scored with four
models: `rf`, `scrambled` (same forest, shuffled training labels), `mean`
(predict the training mean), `size_only` (forest on heavy-atom count alone).

**Why.** The surrogate is the instrument the GA will exploit; making it look
good would hide the phenomenon under study. The controls say what the numbers
mean: scrambled is the H4 control arm's engine and must score ≈ 0;
`size_only` checks the classic confound that bigger molecules bind tighter.

**Result (median [IQR] over 10 seeds).**

| split | RF R² | RF Spearman | scrambled R² | size_only R² |
|---|---|---|---|---|
| random | 0.645 [0.619, 0.659] | 0.811 | −0.122 | 0.050 |
| scaffold | 0.494 [0.477, 0.508] | 0.708 | −0.127 | 0.007 |
| paper | 0.300 [0.255, 0.320] | 0.541 | −0.080 | −0.032 |

RF RMSE: 0.68 / 0.84 / 0.96 log units (random/scaffold/paper), against 1.16–1.21
for the mean predictor. Pooled over seeds, |error| falls as max-Tanimoto to
train rises in every split (Spearman −0.33 / −0.30 / −0.20).

**Interpretation.** Controls behave: scrambled ≈ 0 (negative R² because its
predictions add noise around the mean), so the pipeline is not leaking labels.
Size alone explains ≤5% of variance, so size is not the main driver. Accuracy
degrades cleanly random > scaffold > paper. R² is computed against the
test-set mean.

**Cost.** RF cannot extrapolate beyond its training labels' range and was not
tuned; a mediocre surrogate is fine. Scrambled random-split Spearman is +0.08
[−0.02, 0.11]: consistent with noise at n=403 (SE ≈ 0.05) but worth keeping an eye on.

**How I'd defend it.** "I deliberately did not tune it. I report what it
can do against what random labels and molecule size alone can do, on three
splits and ten seeds."

**What would change my mind.** A scrambled-label Spearman clearly above ~0.15,
or size-only beating the scrambled model by much: either would mean leakage or
a size artifact to chase before continuing.

---

## D-18: One fixed real surrogate and one fixed scrambled surrogate; matched GA runs

**Date:** 2026-09-30
**Stage:** 5 (GA) — **[judgment]**

**Decision.** One random forest trained once on all 2,016 curated molecules
(`random_state = config.SURROGATE_SEED = 42`), and one trained on labels shuffled
once (`SCRAMBLE_SEED = 20260930`), reused unchanged by every GA seed and arm.
GA seeds (42–46) control only the search: start molecules, parent choice and
edits. For each seed, all four arms start from the same 100 training molecules
(drawn at random from those inside the size window) and use identical operators.

**Why.** Separates GA stochasticity from surrogate stochasticity, and makes arms
differ in exactly one thing, the fitness function (paired comparison).

**Alternatives considered.** Retrain the forests per seed (what the first, pilot
run did) — rejected: seed-to-seed variance then mixes search noise with
instrument noise. A surrogate trained on a scaffold-split training set — kept
as a later secondary robustness run, one predefined split.

**Cost.** Conclusions rest on one particular scrambled permutation and one
forest; the scrambled result in particular could differ for another permutation.
Start molecules are training molecules, so generation 0 has max Tanimoto 1.0 and
the forest's predictions for them are in-sample (near their measured values).
Every trajectory therefore begins with an artificial discontinuity at generation 1,
and "best predicted activity" falls below its generation-0 value in three arms
(it is not a failure of optimisation: the best start molecules are known actives
that get displaced).

**How I'd defend it.** "I held the instrument fixed so any difference between
arms is the objective, not the surrogate's random state. Same seed means same
start molecules, so comparisons are paired."

**What would change my mind.** Rerunning with a second scrambled permutation
changing the scrambled-arm conclusion.

---

## D-19: Primary heavy-atom window 20–39; 15–50 kept as a labelled pilot

**Date:** 2026-09-30
**Stage:** 5 (GA) — **[judgment]**

**Decision.** Generated molecules must have 20–39 heavy atoms, the 5th–95th
percentile of the training set (full range 5–78). 15–50 is available
(`--relaxed`) for stress tests only. The first full run, made with 15–50 and
per-seed retrained forests, is preserved unmodified as
`results/05_pilot_relaxed_15-50_populations.csv` (+ console log, + figures
suffixed `_pilot_relaxed`). It is a pilot, not a Stage 5 result.

**Why.** Outside the training size range the surrogate has no examples; any gain
there is extrapolation by construction. Without a window, `activity_only` grew
molecules to ~45 heavy atoms (pilot), beyond anything the surrogate was trained on.

**Cost.** The window truncates the very behaviour (size exploitation) that the
pilot showed. Populations pile up at the bounds: 28% of multi_real and 19% of
druglike_only final molecules sit at 20 or 39 atoms (almost all at 20), so the
drug-likeness objectives are partly bounded by the window, not optimised.

**How I'd defend it.** "I matched the allowed chemistry to what the surrogate
was trained on, so a claim that optimisation left the training distribution
can't be just a size effect. I kept the wider run, labelled, to show what
happens without the constraint."

**What would change my mind.** Bound-piling changing conclusions; then the
window or objectives would need revisiting (e.g. a size-penalty term).

---

## D-20: Four arms, geometric-mean fitness; how Stage 5 is read

**Date:** 2026-09-30
**Stage:** 5 (GA) — **[judgment]**

**Decision.** Arms: `multi_real` (activity × QED × SA, real surrogate),
`multi_scrambled` (same, activity from the scrambled surrogate), `activity_only`,
`druglike_only` (QED × SA, no surrogate). Each term scaled 0–1 (activity linearly
over pActivity 4–9, fixed globally; SA as (10 − SA)/9), combined by geometric
mean. GB-GA: population 100, 100 children/generation, 50 generations, elitist
(keep best 100 of parents + children), crossover + 50% mutation.

**Interpretation rule.** Falling similarity to the training set is **not**
evidence of surrogate exploitation: it falls from 1.0 to ~0.35 in druglike_only,
an arm that never sees the surrogate. Evidence is the movement of activity-driven
arms **relative to druglike_only, same seed** (paired).

**Result (5 seeds, generation 50, median [min, max] over seeds; paired contrast =
arm minus druglike_only).**

| arm | real pred. change vs start | paired Δ real pred. | paired Δ max Tanimoto | median QED / SA / atoms |
|---|---|---|---|---|
| multi_real | +1.21 | +1.85 [+1.47, +1.98] (5/5) | +0.035 [0.00, +0.24] | 0.91 / 2.2 / 21 |
| multi_scrambled | −0.61 | −0.13 [−0.32, +0.44] (1/5) | −0.008 [−0.01, +0.08] | 0.92 / 2.2 / 23 |
| activity_only | +1.91 | +2.50 [+2.37, +2.55] (5/5) | **+0.33** [+0.33, +0.38] | 0.31 / 3.2 / 35 |
| druglike_only | −0.58 | (baseline) | (baseline: 1.0 → 0.355) | 0.94 / 1.7 / 22 |

- Scrambled arm's *own* surrogate rises (+0.5, to ~7.16) while the real surrogate
  does not improve: the optimiser exploits noise but the gain is not transferred.
  H4 holds qualitatively; amplitude of the own-surrogate rise is 0.5 vs 1.2–1.9.
- Activity-driven arms end **closer to** the training set than the baseline, not
  further (paired Δ Tanimoto ≥ 0). Stage 4 calibration (pooled 3 splits, 10 seeds):
  RF MAE 0.39 at similarity ≥ 0.8, 0.51 at 0.7–0.8, 0.67 at 0.6–0.7, 0.75 at 0.4–0.6,
  ~0.9 below 0.4. activity_only molecules mostly sit at 0.6–0.8 (MAE 0.5–0.7),
  i.e. where the surrogate is reasonably reliable; multi_real has ~42% below 0.4.
- Diversity collapse: final populations hold 22/100 (multi_real) and 21/100
  (activity_only) unique scaffolds vs 62–70 for the other arms; activity_only's
  top molecules resemble dinaciclib (the 4KD1 ligand, a high-pActivity training
  compound): local exploitation of known actives, not extrapolation.

**Cost / caveats.** Max Tanimoto is size-dependent (Spearman 0.54 with heavy
atoms, final generation), and activity_only has no molecules at 20–23 atoms, so
its similarity cannot be size-matched from this run. Size-matched at 20–23 atoms,
multi_real 0.41 vs druglike 0.35 vs scrambled 0.36, same conclusion. Stage 4
calibrates on real molecules; GA molecules are naive graph edits, so it may be
optimistic for them. Only 9 Stage-4 test molecules were both predicted ≥ 8 and below
0.6 similarity: essentially no validation data where the GA would exploit.

**How I'd defend it.** "Four arms separate the effect of the objective from the
GA's own drift. The controls behave as designed. The surprise is that activity
pressure keeps molecules near known actives rather than pushing them out of
domain, so I report that instead of the hypothesis I expected."

**What would change my mind.** A size-matched or bit-count-normalised similarity
reversing the activity_only result; or a second scrambled permutation/forest
seed changing H4.

---

## D-21: Headline statistics describe generated molecules only (birth_generation > 0)

**Date:** 2026-09-30
**Stage:** 5 — **[judgment]**

**Decision.** Generation 0 (the 100 starting training molecules) stays in the
trajectory figures as the explicit starting reference. From generation 1 on,
every statistic is computed over population members with `birth_generation > 0`.
Final-generation summaries and counts of "unique generated candidates" use the
same rule. The GA populations and starting molecules are unchanged; the GA was
not rerun.

**Why.** Starting molecules are training molecules: their similarity is 1.0 by
construction and the forest's predictions for them are in-sample. Letting them
into a statistic about generated chemistry mixes the surrogate's memory with its
generalisation. (It is also why "best predicted activity" used to dip below its
generation-0 value: the best start molecules are known actives that get displaced.)

**Cost.** Effect on the final generation is small (starting molecules are 0–2%
of it), but early-generation statistics rest on fewer molecules. Across all
generations 4.8% of population rows are surviving starting molecules.
Superseded numbers: D-20's table included them. Corrected generation-50 medians
(5 seeds): real prediction multi_real 7.75, multi_scrambled 5.91, activity_only
8.60, druglike_only 6.10; best real prediction 8.28 / 6.97 / 8.67 / 6.81; all
paired contrasts vs druglike_only unchanged in sign and within rounding.

**How I'd defend it.** "Statistics about what the GA generated shouldn't include
molecules it was handed."

**What would change my mind.** A need to describe the population as the GA sees
it (survivors included), which `05_ga_populations.csv` still supports.

---

## D-22: Raw ECFP4 Tanimoto stays primary; add a size-conditioned percentile

**Date:** 2026-09-30
**Stage:** 5b — **[judgment]**

**Decision.** Primary similarity = highest Tanimoto between a molecule's ECFP4
fingerprint (Morgan radius 2, 2048 bits, binary, no chirality) and the 2,016
training fingerprints. Second measure, `size_pctile`: that value as a percentile
within the leave-one-out nearest-neighbour Tanimoto distribution of training
molecules whose heavy-atom count is within ±2 of the generated molecule's
(widening by 1 atom if fewer than 100 reference molecules; the rule never
triggered in the 20–39 window: ≥153 per stratum). Both are stored in
`results/05_ga_populations_scored.csv`; the reference is
`results/05_training_loo_similarity.csv`. No corrected-Tanimoto formula is used.

**Why.** Raw similarity depends on size: the training set's own leave-one-out
median is 0.73 for 20–24 atoms and 0.83 for 35–39. A 21-atom molecule is a worse
match than a 35-atom one regardless of novelty. Leave-one-out avoids the trivial
1.0 of comparing a training molecule with itself. 166 training molecules still
have an identical fingerprint to another training molecule (stereoisomers,
possible activity cliffs that ECFP4 cannot distinguish).

**Result (generation 50, generated only, median over seeds).** Size percentile:
multi_real 7.5, multi_scrambled 4.0, activity_only 11.6, druglike_only 4.6.
Every arm ends far below a typical training molecule of its size (50th),
including activity_only, whose raw 0.69 looks close; ordering is unchanged
(activity-driven arms a little closer than baseline: paired +2.1 and +7.5 points).

**Cost.** The percentile is relative to the training set's own analogue-series
density, which is high, so it is a demanding bar. The generation-0 percentile
(~100) is trivial (self-match) and is not a meaningful reference.

**How I'd defend it.** "I kept the plain Tanimoto and added an empirical,
assumption-free way of asking whether it is high or low for a molecule of that
size."

**What would change my mind.** Stage 7: if neither measure separates
well-supported high-prediction candidates from low-similarity candidates with
larger Stage 4 error, the audit should use a different distance.

---

## D-23: How Stage 5 is stated (descriptive, not hypothesis-preserving)

**Date:** 2026-09-30
**Stage:** 5 — **[interpretation]**

**Statement.** (i) druglike_only causes substantial chemical drift by itself
(max Tanimoto 1.0 → 0.36, size percentile ~5). (ii) multi_real raises
real-surrogate prediction strongly (+1.85 over baseline, 5/5 seeds) with about the
same or slightly higher raw similarity than druglike_only. (iii) activity_only
converges on relatively high-similarity potent neighbourhoods while sacrificing
QED and SA. (iv) multi_scrambled optimises its random surrogate (own prediction
+0.5) without transferring that gain to the real surrogate (1/5 seeds above
baseline). Activity optimisation here appears to exploit known potent chemical
neighbourhoods, not to force systematic extrapolation. H2 as worded (similarity
falls as activity climbs) is **not supported** by this experiment.

**Descriptive support (no rerun; `05d_ga_descriptive_checks.py`).**
- Diversity collapse: distinct scaffolds per 100 fall to 22 (multi_real) and 21
  (activity_only) by gen 50, vs 70 and 62. About 48% of activity_only's final
  molecules share one scaffold, 2-anilino-4-(imidazol-5-yl)pyrimidine (31
  training molecules, mean pActivity 8.01); multi_real's top scaffold (32%) is
  the bare imidazolyl-pyrimidine core, which no training molecule has.
  Nearest training neighbours are top-decile actives for 97% (activity_only) and
  62% (multi_real) of final molecules, versus 3–5% in the other arms, and
  activity_only's 492 molecules map to only 14 distinct nearest training molecules.
- Boundary piling is the window capping a real preference: training-set median
  QED falls from 0.74 (20 atoms) to 0.38 (36–39) and SA rises 2.5 → 3.3, so the
  drug-likeness terms push to the lower limit (28% of multi_real final molecules
  at 20 atoms; 58% at ≤21), while within multi_real size and real prediction
  correlate +0.39 and size and QED −0.36. activity_only piles at the upper limit (14% at 39).

**How I'd defend it.** "The surrogate's optimiser learned which known series are
potent and moved toward them, which is what a surrogate trained on those series
should allow. That is a useful negative result for the extrapolation hypothesis."

**What would change my mind.** A stress run with the relaxed window or a
second scrambled permutation giving a different picture.

---

## D-24: Pareto analysis design (Stage 6)

**Date:** 2026-09-30
**Stage:** 6 — **[judgment]**

**Decision.** Objectives: predicted pActivity from the real surrogate (maximise),
QED (maximise), SA score (minimise). Population: unique SMILES among generated
molecules (birth_generation > 0), pooled over all generations and five seeds.
Per-arm fronts, plus one pooled front over all four arms (a molecule found by
several arms counts once, all arms credited). The four molecules drawn are chosen
by a rule fixed before looking: the front's highest prediction, highest QED,
lowest SA, and best GA fitness.

**Why.** The pooled front answers "what trade-offs did the whole search find";
per-arm fronts show which objective produced which part. Pooling over
generations (not just the last one) keeps good molecules that later generations replaced.

**Overlays are stand-ins, and labelled as such.** "Known actives" = top decile
(measured pActivity >= 8.05, 204 molecules) of the TRAINING set: in-sample, not
the held-out positive control (separate experiment, not built). "Random floor" =
300 random training molecules, because no random-ChEMBL sample has been
downloaded; the training set is chemically narrow (kinase-inhibitor-like), so it
is a biased floor. The reference activity axis is measured, the generated one
predicted, so dominance comparisons between them are indicative only.

**Result.** 11,931 distinct generated molecules; pooled front 101 (multi_real 65,
activity_only 23, druglike_only 19, multi_scrambled 5). Front spans predicted
pActivity 5.84-8.68, QED 0.29-0.95, SA 1.38-3.80, 20-37 heavy atoms. No
generated molecule is predicted above 8.89 (the forest's ceiling in practice; the
most potent known actives are measured up to 9.5). 136/204 known actives are
dominated by a front molecule, almost entirely through QED/SA; 20/101 front
molecules are dominated by a known active. Against Stage 4's calibration: 44% of the
front has raw similarity < 0.4 (Stage 4 MAE ~0.91), 26% at 0.4-0.6 (0.75), 21% at
0.6-0.8 (0.57), 10% >= 0.8 (0.39). 40 front molecules are predicted >= 8; 9 of those
have raw similarity < 0.6, the region where Stage 4 has almost no validation data.

**What the drawn molecules show.** The top-predicted molecule is the imidazolyl-
anilinopyrimidine core decorated with a chemically implausible appendage
(N-CH2-S-N linking an enamine): the surrogate keeps predicting high activity
for the core and does not penalise nonsense attached to it. The highest-QED
molecule contains an unusual dihydropyrazine. The best-balanced molecule (20
atoms, 4-isobutylbiphenyl sulfonamide, predicted 7.81 at size percentile 5) is
plausible-looking but far from training chemistry.

**Cost.** No structural-alert or stability filter was applied, so some front
members are not sensible molecules; QED/SA reward resemblance to existing drugs
(partly circular) and cannot see the failure above.

**How I'd defend it.** "The front shows the trade-off the search actually found.
The high-activity end is decorated known chemotypes, the drug-like end is small
generic molecules, and the two ends are separated, which is what you expect if
activity and drug-likeness pull in opposite directions through molecular size."

**What would change my mind.** A held-out known-active set or an external
random-ChEMBL sample changing which comparisons look favourable to generated molecules.

---

## D-25: Applicability-domain audit and structural-alert audit (Stage 7)

**Date:** 2026-09-30
**Stage:** 7 — **[judgment; all thresholds arbitrary and stated]**

**Decision.** Reliability is estimated only from Stage 4's out-of-sample
predictions (12,090 test predictions; 2,014 distinct molecules) against two
measures: raw ECFP4 (Morgan r=2, 2048-bit, binary) max-Tanimoto to the training
set, and the size-conditioned percentile (D-22), whose reference for each Stage 4
test molecule is its own split's training set, leave-one-out. Regions: raw
<0.4 / 0.4-0.6 / 0.6-0.8 / >=0.8; percentile <10 / 10-50 / >=50; predicted
pActivity <6 / 6-7 / 7-8 / >=8. A validation cell is used only if it holds
>= 30 DISTINCT molecules (rows overstate evidence: a molecule recurs across
seeds). "Large error" = |error| >= 1.0 log unit (10-fold). Generated molecules
are classed A (predicted >= 8, raw >= 0.6), B (predicted >= 8, raw < 0.6) or C
(predicted < 8). Structural alerts (RDKit PAINS and Brenk) are recorded, never
used to remove molecules, and kept separate from graph validity and domain support.

**Result 1: error rises as similarity falls, on every split, on one curve.**
Pooled MAE by raw region: <0.4: 0.91; 0.4-0.6: 0.75; 0.6-0.8: 0.57; >=0.8: 0.39
(RMSE 1.10 / 0.93 / 0.74 / 0.52). Random, scaffold and paper splits fall on
essentially the same curve: error is a function of similarity, not of how the
split was made. Spearman(|error|, similarity): raw -0.35, percentile -0.33 pooled;
AUC for flagging |error| >= 1: raw 0.71, percentile 0.69 (paper split 0.61/0.59: weak).
The size percentile adds little beyond raw similarity: within a raw region MAE
barely changes across percentile regions (0.6-0.8: 0.61 / 0.58 / 0.47; 0.4-0.6:
0.75 / 0.75), and the two measures are strongly collinear in the data (most
raw<0.4 molecules are also percentile<10), so their separate value cannot be
tested well here. Raw similarity stays the primary measure.

**Result 2: where Stage 4 can and cannot vouch (distinct validation molecules).**
Predicted >= 8 is validated only at raw >= 0.6 (n=94 at 0.6-0.8, MAE 0.44;
n=54 at >=0.8, MAE 0.36). At raw 0.4-0.6 only 7 distinct molecules are predicted
>= 8 (no estimate) and at raw < 0.4 there are none. At raw < 0.4 and predicted
7-8 (n=89) the forest OVER-predicts by +0.36 on average (MAE 0.95).
**No reliability is claimed for any generated molecule with predicted >= 8 and
raw similarity < 0.6.**

**Result 3: generated molecules on the domain.** Generation-50 share with raw
similarity < 0.4: druglike_only 77%, multi_scrambled 69%, multi_real 43%,
activity_only 0%. **Drift alone is not evidence of exploitation:** druglike_only
never sees the surrogate and drifts furthest. Candidate classes (distinct
molecules, all generations): activity_only 2,059 A / 455 B; multi_real 117 A /
164 B; multi_scrambled 2 A / 1 B; druglike_only 7 A / 0 B. Pareto front (101):
31 A, 9 B, 61 C; 92/101 sit in a cell Stage 4 can estimate. The 9 class B front
molecules are two families: five truncated imidazolyl-pyrimidine analogues at raw
0.53-0.58 (size percentile 7-15, just under the 0.6 boundary, so the count
depends on the threshold) and four diarylalkyl sulfonamides at raw 0.33-0.36
(percentile 3-7), predicted 8.00-8.01, with no training support.

**Result 4: structural alerts (recorded, not filtered).** Brenk: training set
30%, known actives 32%, generated overall 27%; activity_only 56%, multi_real 20%,
druglike_only 12%, multi_scrambled 13%. PAINS 5-11%. High-prediction molecules
carry Brenk alerts at ~54% in BOTH class A and class B (odds ratio 1.01): no
enrichment of alerts with weak support. Alerts are driven by oxygen-nitrogen single
bonds, het-C-het, quaternary N, N-oxide (present in known actives), diketo and thiol.
The Pareto front is nearly alert-free (7%) because QED already penalises
alerts and the front is selected on it; the implausible top-prediction front
molecule is caught by two Brenk alerts (het-C-het_not_in_ring, S-N single bond),
but the 9 class B front molecules mostly carry none. Validity, alerts and domain
support are three different properties and none implies another.

**Cost / caveats.** Stage 4 forests were trained on 80% of the data; the GA
surrogate on 100%, so the calibration is slightly pessimistic for it. Stage 4 tests
real molecules; GA molecules are naive edits, so the curve may be optimistic
for them. The one relevant over-prediction signal (+0.36) rests on 89 molecules
and is pooled across splits. Thresholds (0.6, 8, 30, 1.0) are conventions, not
findings.

**How I'd defend it.** "I measured how wrong the surrogate is as a function of
similarity on molecules it had not seen, then placed the generated molecules on that
curve and refused to make any claim where there was no validation data. Two
measures of distance gave nearly the same answer, so I kept the simpler."

**What would change my mind.** The held-out-actives positive control showing
the surrogate recovers potency in weakly supported regions (then B-class
predictions deserve more trust than Stage 4 suggests), or a different
distance (e.g. forest tree variance) separating large-error molecules much
better than similarity.

---

## D-26: Held-out known-actives positive controls (molecule level and scaffold level)

**Date:** 2026-10-02
**Stage:** 5e-5g — **[judgment; rules fixed before the control GAs were run]**

**Decision.** Two controls, same GA protocol as Stage 5 (4 arms, 5 seeds, 50
generations, 20-39 heavy atoms, same operators). The surrogate is retrained
without the held-out molecules; start molecules come from the remainder; held-out
molecules are used only afterwards, to evaluate.
- **Molecule level:** hold out the 204 top-decile actives (measured pActivity >=
  8.05, the 90th percentile; deterministic, no sampling). Siblings and
  scaffold-mates stay in training. Tests interpolation.
- **Scaffold level:** hold out every scaffold family with >= 5 top-decile
  members, ALL members (potent or not): 9 families, 260 molecules (12.9% of
  the data), 93 of the 204 actives. The other 111 actives stay in training. Tests
  generalisation to unseen ring systems.
Threshold 5 was fixed from family sizes alone (top decile = 100 scaffolds, 76
singletons; thresholds 2-10 all include the two largest families); I had
already seen the main GA converge on one of these families, so the rule is not
blind to Stage 5, but the result does not hinge on the threshold.

**Success criteria (never the surrogate's own prediction).** ECFP4 (Morgan r=2,
2048-bit) Tanimoto to the nearest held-out active, primary threshold 0.6 (0.5-0.8
reported); exact fingerprint matches; whether a generated molecule's nearest
neighbour among all 2,016 is a held-out active, against the chance rate; measured
activity of that neighbour (label transfer, no surrogate); recall of held-out
actives; diversity. Reference rows: the starting molecules and the remaining
training molecules (what siblings alone give).

**Surrogate on the held-out actives, before any GA.** Molecule control: trained
labels reach only 8.05, so predictions cannot exceed it (0% predicted >= 8, median
7.33 vs measured 8.40, bias -1.26) but it ranks them (AUC 0.80 against
out-of-fold predictions of the rest). Scaffold control: median predicted 7.52 vs
8.40, bias -0.94, 11% predicted >= 8, AUC 0.87. Both clearly under-predict the
actives they never saw.

**Support after the holdout (nearest remaining-training Tanimoto of held-out
actives).** Molecule: median 0.72, 88% >= 0.6. Scaffold: median 0.67, 76% >= 0.6,
3% >= 0.8. The scaffold holdout is SOFT: related cores with different scaffold
strings stay in training (e.g. other anilinopyrimidines), so most held-out actives
still have a close relative.

**Results (final generation, generated molecules only, 5 seeds pooled).**

| control | arm | median sim to held-out actives | share >= 0.6 | nearest neighbour is a held-out active (chance) | actives hit (any generation) |
|---|---|---|---|---|---|
| molecule | activity_only | 0.71 | 83% | 70% (10%) | 112/204 |
| molecule | multi_real | 0.32 | 1% | 24% | 72/204 |
| molecule | multi_scrambled / druglike_only | 0.28 / 0.29 | 0% / 0% | 8% / 9% | 30 / 16 |
| molecule | remaining training (no search) | 0.33 | 26% | 7% | ceiling: 88% have a sibling >= 0.6 |
| scaffold | activity_only | 0.45 | 22% | 30% (5%) | 56/93 |
| scaffold | multi_real | 0.30 | 1% | 12% | 30/93 |
| scaffold | multi_scrambled / druglike_only | 0.17 / 0.22 | 0% / 0% | 1% / 2% | 22 / 15 |
| scaffold | remaining training (no search) | 0.21 | 6% | 2% | ceiling: 76% |

- Exact rediscovery is rare (<= 0.2% of molecules at identical fingerprint).
- Activity regime by measured label transfer: among generated molecules with a
  neighbour >= 0.6, 75% (molecule) and 96% (scaffold) of activity_only's neighbours
  are themselves top-decile actives.
- Diversity: activity_only's 408 (molecule) / 110 (scaffold) near-held-out
  molecules map to only 10 / 3 distinct held-out actives and 33 / 35 scaffolds;
  top scaffold share 71% / 44%. The same collapse as Stage 5.
- Scaffold families (9): exact scaffold generated by any arm for 6; but 3 of those are
  tiny generic scaffolds (13-18 atoms) that every arm including druglike_only
  produces. Non-generic families found exactly: activity_only 2 (a 20-atom
  thiazolo-indole imine and a 26-atom benzimidazolyl-pyrimidine), multi_real 1
  (a 22-atom pyrazolo-triazine); 3 of 9 never. The raw count "5 of 9" overstates.

**Interpretation.** (1) activity_only passes the molecule-level control clearly (83%
within 0.6, 7x the chance nearest-neighbour rate, far above the no-search baseline
of 26%), but that is sibling-supported interpolation: 88% of held-out actives had
a training sibling to begin with, and recall (55%) stays below that ceiling.
(2) At scaffold level activity_only still finds chemistry near the held-out potent
families (22% >= 0.6 vs 6% baseline; 6x chance; 60% recall vs 76% ceiling; 2 non-
generic scaffolds recovered exactly), which is partial support for generalisation, but
because the holdout is soft it is not evidence of broad generalisation to genuinely
unseen chemistry. (3) multi_real, the main multi-objective arm, FAILS both
controls: it is at or below the no-search baseline (1% >= 0.6) because QED and SA
pull molecules small and drug-like and away from large potent ones. This is
informative, not a pipeline bug: the Stage 5 "exploits known potent neighbourhoods"
reading holds for the activity-only objective, not for the composite. (4) The
negative arms behave as negative controls (0% >= 0.6, well below baseline).

**How I'd defend it.** "I removed known potent molecules, then removed whole potent
series, retrained, and checked whether the search found chemistry near them,
judged against measured activity, with no-search and scrambled/drug-likeness
baselines. The activity-driven search finds them; the composite objective does
not; and I say the scaffold holdout is soft."

**What would change my mind.** A stricter scaffold holdout (e.g. also removing
molecules above a similarity threshold to the held-out series) changing the picture;
or other split seeds / hold-out thresholds giving a different ordering of arms.

---

## D-27: Secondary robustness run: surrogate trained on a scaffold split's training set

**Date:** 2026-10-02
**Stage:** 5h — **[judgment]**

**Decision.** One predefined run (Stage 3 scaffold split, seed 42; 1,613 training
molecules, 403 held out), same GA protocol as Stage 5 (4 arms, 5 GA seeds, 50
generations, 20-39 heavy atoms), with the real and scrambled forests trained on
those 1,613 molecules only. Similarity is to those 1,613. Compared with the
primary run on six claims fixed in `05h_scaffold_surrogate_robustness.py`
before running it. More split seeds only if the picture proved unstable.

**Why.** The primary surrogate saw every molecule, so "training distribution" and
"all known chemistry" coincide. The scaffold-trained surrogate has a real
unseen set (accuracy there: R2 0.51, RMSE 0.84, Spearman 0.70), a check that the Stage 5
picture is not an artifact of training on everything.

**Result (generation 50, generated molecules, median over 5 seeds).** Five of the six
pre-written claims hold under both surrogates:
- multi_real raises the real prediction above druglike_only in 5/5 seeds (both).
- activity_only ends closer to its training set than druglike_only in 5/5 seeds (both;
  similarity 0.79 vs 0.35 in the secondary run, with QED 0.23 and SA 3.9).
- multi_scrambled beats druglike_only on the real prediction in 1/5 (primary) and 0/5 seeds.
- druglike_only drifts as far as anything: 78% of its molecules below raw similarity 0.4,
  vs 0% for activity_only (both).
- multi_real and activity_only keep fewer distinct scaffolds than druglike_only
  (22/21 vs 62 primary; 25/36 vs 40 secondary).
**One claim fails by the pre-written rule:** "multi_real raw similarity >= druglike_only"
holds in the primary run (+0.03) and fails in the secondary (-0.01).

**Honest reading.** The failure is by 0.01 on a median of 0.34 vs 0.35 with tight
ranges, so the robust statement is "multi_real's similarity is about the same as
druglike_only's", which is how D-23 phrased it; "slightly higher" is not robust. In
the primary run multi_real's similarity varied a lot between seeds (share below 0.4
ranged 0% to 88%), in the secondary it was uniformly low (82-91%). Two soft spots:
the diversity collapse is robust for multi_real but marginal for activity_only in
the secondary run (36 vs 40 scaffolds per 100, ranges overlapping heavily), and
druglike_only's own diversity differs between surrogates (62 vs 40). I did not
change the rule to make the claim pass.

**Cost.** One split seed only; claims about arm differences rest on 5 GA seeds, not on
independent splits. Absolute similarity and prediction values are not comparable across
the two experiments (different training sets, different surrogates).

**How I'd defend it.** "The main conclusions (activity pressure raises predicted activity,
activity-driven search ends close to known chemistry rather than far from it, the scrambled
surrogate's gain doesn't transfer, drug-likeness alone drifts as far as anything) hold
when the surrogate has truly unseen test chemistry. One detail, whether multi_real is
slightly more similar than the baseline, does not, and I state it as 'about the same'."

**What would change my mind.** More split seeds (43, 44, ...) flipping an arm ordering, or
the multi_real similarity difference growing in either direction.

---

## D-28: Robustness run extended to five scaffold splits (supersedes the single-split reading in D-27)

**Date:** 2026-10-02
**Stage:** 5h — **[judgment; stability rule fixed before the extra splits were run]**

**Decision.** The scaffold-trained GA was repeated for split seeds 42-46 (the first five
predetermined split seeds), each with the same 5 GA seeds: 25 (split, GA seed) pairs
besides the primary run. Rule declared in advance: a claim is STABLE if it holds in the
primary run and in at least 4 of the 5 splits, otherwise FRAGILE. The six claims are those of D-27.

**Result.**
| claim | primary | splits holding | verdict |
|---|---|---|---|
| multi_real raises the real prediction above druglike_only in every GA seed | holds | 5/5 | STABLE |
| activity_only ends closer to training than druglike_only in every GA seed | holds | 5/5 | STABLE |
| druglike_only drifts at least as far as activity_only (share < 0.4: 84% vs 0%) | holds | 5/5 | STABLE |
| multi_real raw similarity >= druglike_only (median paired difference) | holds | 3/5 | FRAGILE |
| multi_scrambled beats druglike_only on the real prediction in at most 2 GA seeds | holds | 3/5 | FRAGILE |
| multi_real and activity_only keep fewer distinct scaffolds than druglike_only | holds | 3/5 | FRAGILE |

**What the fragile claims actually say.**
- *multi_real vs druglike_only similarity:* per-split paired difference ranges -0.015 to
  +0.048; pooled over 25 pairs median +0.014, IQR [-0.014, +0.041], positive in 15/25. The
  robust statement is that multi_real's similarity is **about the same** as the
  baseline's. "Slightly higher" is not supported.
- *multi_scrambled on the real prediction:* the rule asked for <= 2 of 5 seeds above
  the baseline. Under a true null effect about half the seeds would exceed it, so this
  rule was too strict for a null; the claim as written is fragile (splits 44 and 46 break it,
  scrambled above baseline in 5/5 and 3/5 seeds). Judged by effect size the picture is
  clear: pooled median difference **-0.15** (IQR [-0.23, +0.07], above baseline in 10/25),
  against **+1.55** (minimum +0.86) for multi_real over the same baseline. The scrambled
  arm sits at the baseline on the real surrogate (one outlier of +1.14 in split 46); it
  does not transfer a real gain. The declared rule is left as written and the
  effect-size reading is added beside it, not instead of it.
- *diversity:* multi_real has fewer scaffolds than druglike_only in all five splits
  (23-33 vs 34-50 per 100), a robust collapse. activity_only's collapse is fragile
  (pooled median -11 scaffolds per 100, IQR [-21, +10], fewer in 16/25 pairs; in splits
  45 and 46 it is as diverse as or more diverse than the baseline). The primary run's
  "activity_only collapses to 21 scaffolds" does not generalise across splits.

**Cost.** Five splits share overlapping molecules, so they are not independent draws; the
arm comparisons rest on 5 GA seeds per split. Absolute similarity and prediction values
are not comparable between the primary and scaffold-trained experiments.

**How I'd defend it.** "I declared the rule in advance and reported it as written: three of
six claims are stable, three fragile. The headline conclusions are the stable ones; for the
fragile ones I state the weaker claim the data supports ('about the same', 'no transfer',
'multi_real collapses, activity_only only sometimes')."

**What would change my mind.** More splits moving the pooled differences clearly away
from their present near-zero values.

---

## D-29: Docking setup and the redocking control (Stage 8a-8b)

**Date:** 2026-10-02
**Stage:** 8 (docking; time-boxed, droppable) — **[mixed: standard tooling, judgment on box and protonation]**

**Decision.** Receptor: PDB 4KD1 chain A only. Ligand 1QK, EDO and all waters removed;
hydrogens added at pH 7.4 with pdbfixer/OpenMM defaults (histidine tautomers left to
OpenMM's heuristic, no manual flips); converted to PDBQT with Meeko (never ADFR or
MGLTools). Search box: a 22 A cube centred on the crystal ligand's atoms. Vina 1.2.7,
exhaustiveness 8 (default), 9 poses. The native ligand is docked as deposited
(1-hydroxypyridinium cation, N-OH, hydrogens added with generated coordinates).
Installed `gemmi` (conda-forge), the missing dependency that kept Meeko from importing.

**Why.** 4KD1 is complete (0 missing residues or atoms), monomeric, unmutated (D-08).
The box is large enough for any ligand we generate (29-heavy-atom crystal ligand spans
about 9 A) yet small enough that Vina's search stays focused on the ATP site.

**Redocking result (control 7).** Top-ranked pose RMSD to the crystal pose, heavy atoms,
symmetry-aware: **0.63 / 0.65 / 0.65 A** for Vina seeds 42 / 43 / 44, against the 2.0 A
limit. Top pose scores -9.48 / -9.44 / -9.48 kcal/mol, essentially identical across seeds
(engine noise, control 9, is small here). The crystal pose scores -8.46 as deposited and
-9.45 after local relaxation, i.e. the redocked pose is as good as the relaxed crystal pose,
so the scoring function does not prefer a wrong pose. About 8 s per docking (8 cores).

**Cost / caveats.** This is the easiest possible test: self-docking into the structure
crystallised with that very ligand, so the pocket already has the right shape (no induced
fit to predict). Passing is necessary, not sufficient; it says nothing about how Vina
treats ligands unlike dinaciclib. Rigid receptor, no crystal waters, no cyclin A (D-02),
heuristic histidine protonation, no tautomer or protonation enumeration for ligands.
Vina scores are not binding free energies, and the score correlates with ligand size
(handled by the planned heavy-atom-count baseline).

**How I'd defend it.** "I first checked that the setup reproduces the crystal binding mode
within 0.7 A over three random seeds. I know that is the easy case and I say so; the docking
results are only used as an orthogonal check with size and engine-noise controls."

**What would change my mind.** A cross-docking test (another known CDK2 ligand into 4KD1, or
dinaciclib into 1KE5) failing badly.

---

## D-30: Docking the molecule sets, and what the Vina scores can and cannot tell us (Stages 8c-8g, 9)

**Date:** 2026-10-02
**Stage:** 8-9 — **[judgment; set sizes chosen by the author, all selection seeded]**

**Decision.** Dock 788 distinct molecules, one Vina run each (seed 42, exhaustiveness 8):
100 known actives (random from the top decile, measured pActivity >= 8.05); 100 random
training molecules; 100 from the final generation of each of the four GA arms; the whole
pooled Pareto front (101); and 100 property-matched decoys (heavy atoms, MW, cLogP, H-bond
donors/acceptors, rotatable bonds; candidates from the druglike_only and multi_scrambled
arms, each dissimilar to every top-decile active at ECFP4 Tanimoto < 0.5). Thirty
molecules were re-docked with two more Vina seeds. This took ~105 minutes, over the
project's ~2-hour compute budget (the author chose the "fuller" option knowingly).
One molecule failed to prepare. Matching quality of decoys: heavy atoms SMD +0.02, but MW
-0.43 and H-bond donors -0.44 (decoys are a bit lighter and have fewer donors).

**Result 1: Vina does not separate known actives from the references.**
AUC (actives vs decoys) 0.53, vs random training 0.45; heavy-atom count alone 0.48 and 0.37
(actives are smaller). Spearman between -Vina and measured pActivity over actives + random
(n=200): +0.02 (size-adjusted +0.03); +0.14 within the random set alone. The redocking
control (D-29) passed, so the setup is not broken, but on this chemistry the score carries
essentially no activity signal. A pose-based feature does a little better: heavy N/O
within 3.5 A of the hinge backbone (Glu81 O, Leu83 N or O) gives AUC 0.58 (vs decoys) and
0.62 (vs random), Spearman +0.19 with measured pActivity within the random set.

**Result 2: size.** Over all 788 molecules, Vina improves by 0.06 kcal/mol per heavy atom
(Spearman -0.49). The size line fitted on the random+decoy sets alone is much shallower
(-0.017) because those sets span a narrow range, so the stratified comparison (median score
within 20-24, 25-29, 30-34, 35-39 heavy atoms) is the safer control; both are in 08e.

**Result 3: H3 (do optimized molecules beat known actives on Vina?).** Probability that a
molecule of the set scores better than a random known active (0.5 = same): activity_only
**0.73** [0.66, 0.80] (size-adjusted 0.70; within size strata 0.61-0.75), multi_real 0.33,
multi_scrambled 0.29, druglike_only 0.26, Pareto front 0.45 [0.37, 0.53]. So the
multi-objective optimized molecules and the front do NOT beat known actives, as H3 says, but
the activity-only molecules DO, in every size stratum where both exist. Because Vina does not
rank known actives (Result 1), "scores better" is not evidence of real potency: **H3 is
inconclusive as an orthogonal test, since the orthogonal signal has no demonstrated
validity on this chemistry.**

**Result 4: does the surrogate's predicted gain show up in docking?** Spearman(predicted
pActivity, -Vina) within a set: multi_real -0.23, multi_scrambled +0.09, activity_only +0.28,
druglike_only +0.03, Pareto front +0.62 (the front spans two ends, drug-like molecules at
low prediction and activity-driven ones at high prediction, so it mixes two populations
and is not a within-chemistry correlation).

**Result 5: hinge contacts and poses.** Share of poses with a hinge contact: activity_only 96%,
known actives 89%, decoys 81%, random training 68%, Pareto front 62%, multi_scrambled 58%,
multi_real 55%, druglike_only 44% (crystal dinaciclib: 2.65 A). Activity pressure yields
hinge-binding chemistry, consistent with D-23; the drug-likeness-only arm loses it.
But 81% of decoys also touch the hinge: decoys are mutated descendants of kinase-inhibitor-like
molecules, not clean negatives, so the decoy test is pessimistic for any method.
PoseBusters: 99% of 787 poses pass every check (all sets 100% except activity_only 94%, six poses
flagged as radicals from hypervalent sulfur in generated molecules, and one decoy). Poses from a
docking engine are physically sane by construction, so a high pass rate is weak evidence; the
failures point to odd generated molecules, not bad poses.

**Result 6: engine noise (control 9).** 30 replicate molecules, 3 seeds: median SD 0.01
kcal/mol, maximum range 1.03; rank agreement between seeds Spearman 0.93 and 0.92.

**Cost / caveats.** Rigid holo receptor, no waters, no cyclin A, heuristic histidine
protonation, single protonation/tautomer per ligand, unspecified stereocentres set by the
embedding, Vina is not a binding free energy. The random reference is a stand-in for the
planned random-ChEMBL sample. One docking per molecule at one seed in the main analysis.

**How I'd defend it.** "I docked known actives and property-matched decoys to ask whether Vina
is a usable orthogonal check at all. It is not, here: it does not separate them and does not
track measured activity, with size and engine noise controlled. So I do not use docking
as evidence for or against the surrogate's gains; what the poses do show is that the
activity-driven search finds hinge-binding chemistry."

**What would change my mind.** A better docking or rescoring setup (cross-docking several CDK2
structures, per-series analysis, consensus scoring) that does separate actives from decoys.

---

## D-31: H5 mitigation arm: design fixed before running (Stage 10)

**Date:** 2026-10-02
**Stage:** 10 — **[judgment; fixed before any constrained run]**

**H5.** Constraining the GA to stay near the training distribution attenuates the activity gain
(H1) but improves transfer. Measure the exchange rate.

**Decision.**
- *Constraint:* a hard similarity floor. A generated molecule whose raw ECFP4 (Morgan r=2, 2048-bit)
  max-Tanimoto to the training set is below tau has fitness 0, so it cannot survive selection.
  tau in {0.4, 0.5, 0.6, 0.7}; tau = 0.6 is the boundary below which Stage 7 found the surrogate's
  high predictions unvalidated. Unconstrained = the existing Stage 5 / 5f runs (no new unconstrained runs).
- *Arms constrained:* multi_real (the composite objective) and activity_only. Everything else is
  as in Stage 5: same GA, seeds 42-46, 50 generations, 20-39 heavy atoms, same start molecules per
  seed, so each constrained run is paired with its unconstrained twin.
- *Where:* (a) the primary surrogate (all 2,016 molecules), to measure attenuation; (b) the molecule-level
  and (c) the scaffold-level held-out-actives controls (D-26), to measure transfer. Their surrogates
  are those of D-26; "training set" in the constraint is the control's remaining training molecules.
- *Cost axis (H1):* median real-surrogate prediction of generated final-generation molecules.
- *Transfer axis:* measured, not surrogate-based: share of final-generation generated molecules within
  Tanimoto 0.6 of a held-out active, share whose nearest neighbour among all 2,016 is a held-out active,
  and recall of held-out actives (any generation, pooled seeds). Docking is NOT used as a transfer
  measure: D-30 found Vina does not discriminate actives from decoys on this chemistry.
- *Exchange rate:* for each arm and tau, change in recovery (percentage points of molecules within 0.6
  of a held-out active) per 1.0 pActivity of predicted activity given up, versus the unconstrained twin.

**Why.** The held-out controls give experimental ground truth for "transfer"; the similarity floor is
the simplest possible constraint and ties to the reliability analysis of Stage 7.

**What would count as support.** Predicted activity falls as tau rises (attenuation) AND recovery of
held-out actives rises. Attenuation without a transfer gain, or no attenuation, is a result
against or beyond H5, and will be reported as such.

**Cost / caveat.** A hard similarity floor to the training set also shrinks the space of novel molecules,
so any "transfer" gain may simply be the constraint keeping molecules near known actives (in the
controls, near the siblings of the held-out ones); the scaffold control is the less trivial test.

---

## D-32: H5 result: a similarity floor does not attenuate predicted activity, and does not reliably improve transfer

**Date:** 2026-10-02
**Stage:** 10 — **[result; design in D-31]**

**Setup as declared.** Hard floor tau in {0.4, 0.5, 0.6, 0.7} on raw ECFP4 max-Tanimoto to the
training set; arms multi_real and activity_only; 3 experiments x 5 seeds x 4 floors x 2 arms = 120
runs, each paired by seed with its unconstrained twin; transfer measured on the two held-out-actives
controls, not by docking (D-30).

**Result A, cost (primary surrogate).** Predicted activity is barely touched. multi_real: 7.75 unconstrained
vs 7.88 / 8.11 / 7.85 / 8.13 at tau 0.4-0.7 (paired change +0.09, +0.07, 0.00, +0.17, never reliably
negative); activity_only: 8.60 vs 8.55 / 8.53 / 8.58 / 8.50 (paired change 0.01 to -0.07). The cost of
the floor appears elsewhere: for multi_real, QED falls from 0.91 to 0.72 at tau 0.7, SA rises from 2.2 to
2.6-2.8, and molecules grow from 21 to 26 heavy atoms; activity_only is unchanged (QED 0.3, 33-35 atoms).
The "attenuates H1" half of H5 is **not supported**: high predictions live near the training data.

**Result B, transfer (held-out controls; median over 5 seeds [min, max]).**
- Molecule level, multi_real: share of generated molecules within 0.6 of a held-out active rises from
  0% to 19 / 48 / 74 / 77% with the floor (+67 pp at tau 0.6); NN-is-held-out 24% to 30-52%; recall 35% to
  40-53%. But the floor keeps molecules near the TRAINING set, which contains siblings of the held-out
  actives (88% have one at >= 0.6), so this gain is largely trivial. activity_only already recovers
  96% [36, 100] unconstrained and gains nothing (69-98%, noisy).
- Scaffold level (the less trivial test): multi_real 1% to 3-6% (ranges up to 23%): no meaningful gain;
  activity_only 28% [0, 52] to 0-10% (it falls, with enormous seed-to-seed spread, 0 to 84%). Recall pooled
  over all generations does rise (multi_real 32% to 66% at tau 0.7, activity_only 60% to 60-68%), reflecting
  early-generation molecules kept near the start rather than better final molecules.

**Result C, exchange rate.** Undefined for 13 of 16 conditions: the floor does not lower predicted
activity (it rises or is flat), so there is nothing "given up". The few defined rates are negative
(recovery and activity both lower, activity_only). The honest exchange is not against activity but against
drug-likeness: for multi_real at the molecule level the floor buys about +67 pp recovery for -0.06 QED
(tau 0.6), at the scaffold level +2 pp for the same -0.06 QED.

**Interpretation.** H5 as worded (floor attenuates the gain, improves transfer) is not supported:
there is no attenuation, and transfer improves only where it can be explained by proximity to siblings
(molecule level, composite arm). The surrogate's high predictions are already in well-supported
chemistry (Stages 5 and 7), so keeping the search there costs little predicted activity but pulls the
composite objective toward larger, less drug-like molecules.

**Cost / caveats.** Five GA seeds with very wide ranges (activity_only scaffold recovery 0-84%); a hard
floor on one similarity definition; transfer judged at one threshold (0.6); the controls' training sets
include siblings of the targets (molecule level) and related cores (scaffold level, 76% of held-out actives
still have a relative at >= 0.6). Docking was not used as a transfer measure because it carries no signal.

**How I'd defend it.** "I declared the constraint and the metrics before running, and measured transfer
against held-out measured actives, not the surrogate. The floor was nearly free in predicted activity, which
tells me the surrogate was already being used where it is supported; where it improved recovery it did so
by staying near siblings of the targets."

**What would change my mind.** A soft penalty (rather than a hard floor), a different similarity measure,
or more seeds narrowing the activity_only scaffold-level range and showing a gain.

---

# Repair pass of 2026-10-03

Entries D-33 to D-39 correct implementation defects and imprecise wording found in review. They do not delete or edit
earlier entries; where an earlier statement is superseded, the entry says which. Historical results files were not overwritten.
No full GA or docking campaign was rerun; the only new computation at scale was the read-only PoseBusters re-audit of the
saved poses. Smoke tests of new code paths are not scientific results.

## D-33: What kind of change each repair is

**Date:** 2026-10-03 — **[policy]**

A. *Implementation repairs that preserve the experiment:* similarity floor as a feasibility filter (D-36), output-path protection and
manifests (D-36), resume-safe docking (D-35). Checked by reproducing saved populations exactly (generations 0-3, seed 42, `05_ga_populations.csv`
and the tau = 0.6 constrained file) under the legacy policy.
B. *Corrected analyses, written beside the historical files:* PoseBusters re-audit (D-34), positive-control evaluation with explicit definitions
(D-38), leakage and artifact audits (D-39), corrected figure 5 and Table 1.
C. *New experiment variants that still need full runs:* the `corrected` candidate-preparation policy with a closed-shell rule (D-37). Only
mechanism smoke tests exist.
None of this establishes biological activity, fixes H3, or validates H5 transfer.

## D-34: PoseBusters aggregation was incomplete; corrected audit (supersedes the pass rates in D-30)

**Date:** 2026-10-03 — **[correction]**

**Defect.** `scripts/08f` chose check columns with `dtype == bool`. A check holding True/False and any missing value is object-typed and was
dropped from the audit for every molecule; `internal_energy` was dropped that way (the saved report has 21 of the 22 dock-configuration checks).
**Cause, established from PoseBusters 0.6.5's own code:** `energy_ratio.py` returns NaN (including for `energy_ratio_passes`) when a molecule has no conformer,
does not sanitize, lacks force-field parameters, cannot be converted to InChI or no ensemble can be built; `Energy ratio module failed` messages were logged during the
first audit. A check that is merely inapplicable is NOT reported as NaN (with no cofactors or waters, those checks return True), so a missing value always means
"not evaluated". The saved artifacts alone could not establish this; the library source and a run on three poses did.
**Repair.** `posebusters_audit.py` freezes the 22 required check names (tested against the installed API), keeps pass, fail and missing apart per check, defines status
(fail = any explicit False; unevaluable = no False but some missing; pass = all True), keeps missing counts for failed poses, identifies poses by stable ids passed
through PoseBusters in the molecule name, reads SDF records by position (a malformed record keeps its place), validates id uniqueness and join cardinality, and counts
molecules without an audit row as not passing.
**Result on the 787 saved poses (`results/posebusters_audit/`):** 774 pass, 13 fail, 0 unevaluable. Failed checks: no_radicals 6, internal_energy 6,
aromatic_ring_flatness 1, bond_angles 1. Five poses have a missing `internal_energy` result and also fail another check (they are counted as failures; their missing count
is retained). Six poses the first audit passed fail `internal_energy`. Pass rates (pass / n in set): known_actives 100%, random_training 100%, decoys 98% (1 molecule failed to dock,
so has no pose), multi_real 98%, multi_scrambled 100%, activity_only 92%, druglike_only 100%, pareto_front 98%.
**Superseded wording.** D-30's "99% of poses pass every check" and "100% for every set except activity_only 94% and one decoy": read 98.3% and the figures above. Historical files
`results/08f_posebusters_by_molecule.csv`, `figures/08_posebusters.png`, `figures/final/fig5_orthogonal_validation.png`, `results/table1.csv` are kept and use the 21-check audit;
corrected variants are `fig5_orthogonal_validation_corrected_pb.png` and `table1_corrected_pb.*`. A high pass rate remains weak evidence: a docking engine places ligands without
overlaps by construction.

## D-35: Resume-safe docking job store

**Date:** 2026-10-03 — **[implementation repair]**

**Defect.** The first `08d` rewrote a pose file when a resumed process first wrote for a seed, trusted existing CSV rows without checking poses, could write `status="ok"` with
no pose, and reused any grid-map directory. **A read-only audit found no damage:** for every Vina seed the saved score rows and pose records agree exactly (787 / 30 / 30
ok scores and poses, no duplicates, no orphans; the one preparation failure is recorded), so the historical docking results are not claimed corrupted.
**Repair (`dock_jobs.py`, `docking.py`, `08d`).** A job is molecule x Vina seed x ligand-preparation seed; its pose and record are separate files written atomically, record last, so a
record exists only for a finished job. `ok` requires a finite score and a readable pose; failures are `prep_failed`, `dock_failed`, `pose_export_failed`. Failures are not retried unless
`--retry-failed`. On resume, orphan poses, missing or unreadable poses, incomplete, mismatched or duplicate records are reported and block the run unless `--repair`. A manifest holds the
receptor hash, box, Vina settings, preparation seed and versions; a mismatch refuses the resume with the differing key; grid maps are accepted only with a hash manifest and every file
present and matching. Legacy artifacts (no manifest) are never adopted. New runs write to `--run-dir`; the legacy `receptor_maps/`, `poses_seed*.sdf` and `08d_docking_scores.csv`
are untouched. A smoke run of three real dockings resumed with zero new dockings.

## D-36: Output safety, manifests, and the similarity floor as a feasibility constraint

**Date:** 2026-10-03 — **[implementation repair]**

`--relaxed` used to keep the primary output name; now the primary name is reserved for the primary protocol (legacy policy, 20-39 window, default arms, seeds, generations),
other runs derive an explicit name (e.g. `05_ga_relaxed_populations.csv`, `05_ga_seeds-42_gen-3_populations.csv`), `--out` is honoured, existing files are never replaced without
`--overwrite`, and every new run writes `<output>.manifest.json` (settings, package versions, platform, git commit, rejection counts). The same applies to `05f` and `10`.
The similarity floor used to set fitness 0 below the threshold, which does not exclude a molecule when too few positive-fitness candidates exist. It is now a filter before survivor
selection: starting molecules must satisfy it (else an error; the floor is never relaxed), invalid or non-finite thresholds and non-finite similarities are handled explicitly,
feasible zero-fitness molecules stay distinct from infeasible ones and survive when needed. The saved constrained populations contained no floor violations (audited read-only);
this is a robustness repair, not evidence against them.

## D-37: Corrected candidate preparation and closed-shell acceptance (a NEW experiment variant)

**Date:** 2026-10-03 — **[new variant; not run at scale]**

Training molecules were cleaned, reduced to the largest fragment, neutralised and tautomer-canonicalized; generated molecules were scored as edited. `--policy legacy` (default, all
historical runs) keeps that. `--policy corrected` standardizes each candidate with the training policy, validates the result (sanitizes, one fragment, size window, no radical
electrons), deduplicates after standardization, scores activity, similarity, QED and SA on that one stored structure, keeps the pre-standardization proposal and the parents for
traceability, and refuses an initial population that standardization would shrink. Rejection reasons are counted per run. The legacy acceptance admits radicals (saved activity-only
populations contain six flagged by PoseBusters); the corrected policy rejects them. This changes the search, so its results are a different experiment and must never be merged
with legacy results or used to update earlier conclusions. Limits: passing these checks does not show stability or synthesizability. **Tautomer enumeration is incomplete for
23 of the 2,016 curated molecules (1.1%, "MaxTransformsReached", default limit 1000), so tautomer invariance is not guaranteed for them**; the status is recorded and
counted ("tautomer_incomplete"). RDKit's tautomer step also strips sp3 and double-bond stereo where tautomerism could interconvert it (amino-acid enantiomers merge), in the
project's curation as well. Two smoke runs (1 seed, 2 generations) showed the mechanism works and the policies are distinguishable (6 of 105 generated structures differ between
proposal and scored form); they are not results.

## D-38: Positive-control evaluation, corrected definitions (supersedes labels in D-26 and D-32)

**Date:** 2026-10-03 — **[correction of labels; conclusions unchanged]**

Problems in the first evaluation: fingerprint equality was called "exact" recovery; a nearest neighbour's measured activity was described as activity transfer; starting molecules
matched themselves (nearest-held-out share 0%) while the remaining-training baseline excluded self-matches; and n held-out / n total was called a chance rate.
Corrected (`recovery.py`, `05g` and `10b` writing `_v2` outputs, historical files kept): *proximity* (ECFP4 Tanimoto >= 0.6 to a held-out active) replaces "recovery"; fingerprint identity is
separate from exact identity, which is reported under three keys (standardized isomeric, standardized flat, strict stereo without tautomer merging); own identity is excluded for starting
molecules and the remaining-training baseline; generated molecules are not excluded; nearest-neighbour ties go to the lowest index and are counted; the old "chance rate" is only the
reference-set prevalence of held-out actives and carries no inferential claim (no replacement null model is proposed; the no-search rows are the empirical comparators); a neighbour's
measured pActivity is an annotation of that neighbour; weighted and distinct counts are both given.
Numbers that changed: starting molecules' nearest-held-out share 8% (molecule control) and 1% (scaffold control), not 0%; exact standardized identity with a held-out active is 0.2% of
activity-only molecules at the molecule level (about one molecule) and 0 elsewhere. The proximity shares, recall and the H5 trade-off values are unchanged apart from names. D-26's
"7x the chance nearest-neighbour rate" and D-32's "recovery" should be read as: 70% vs 7-8% for the no-search baselines, and "proximity share".

## D-39: Leakage and artifact audits; wording corrections

**Date:** 2026-10-03 — **[audit; corrections]**

Read-only audit of the stored splits (`results/audit/`; no split regenerated; alignment of every split to the curated rows verified). Full standardized identity never crosses a split
(duplicates were merged by InChIKey). Stereoisomer-level (connectivity-block) overlap does, test molecules overlapping training, median over 10 seeds [min, max]: random 23 [16, 29],
scaffold 0, paper 3 [1, 10]; identical ECFP4 fingerprint: 27 [21, 33], 0 [0, 1], 3 [1, 10]; exact scaffold: 266, 0, 84; any shared document: 388, 380, 113 (primary document only: 382, 357, 0).
Positive controls: 4 of 204 molecule-level held-out actives have an identical-fingerprint twin in training (stereo or connectivity variants), none at scaffold level; every held-out molecule
in both controls shares at least one document with the remaining training molecules. Docking score-to-pose coverage: complete for every Vina seed. PoseBusters check coverage: 21 of 22 (D-34).
Wording corrected (earlier text stays in history): `primary_document` is the first non-null document met when records of a structure were aggregated, not a verified earliest publication, so
"first paper" in D-13 and D-14 is not verified and the paper split is not a temporal split; "ECFP4 cannot represent activity cliffs" (D-12, CLAUDE.md) should read that molecules with
identical fingerprints (e.g. stereoisomers; 166 curated molecules have an identical-fingerprint twin) cannot be told apart and similarity-based models tend to smooth over cliffs; "Vina has no
entropy" (D-29, D-30) should read that it has an approximate rotatable-bond penalty but no rigorous binding-entropy treatment; the 20-39 heavy-atom window is the central 5th-95th percentile
(the training set spans 5-78), not its full support; elitism preserves molecules by scalar fitness, not every objective component; "newly proposed" means not in the current population or batch,
not globally novel (no novelty filter was added, since that would change the search); "starting molecules had in-sample predictions" still holds. The hinge metric (D-30) is a distance-only
N/O proximity within 3.5 A, not a verified hydrogen bond.

## D-40: Campaign mechanism: a separate directory per run, two preparation policies, nothing silently falls back to historical files

**Date:** 2026-10-04
**Stage:** all (repair and rerun, campaign `c2026-10-03`)

**Decision.** Every stage reads and writes inside `campaigns/<id>/{shared,legacy,corrected}/` chosen by environment variables (`CDK2_CAMPAIGN`, `CDK2_SCOPE`); with no campaign selected
`config.py` points at a sentinel and `require_campaign()` stops the script. Historical results (`data/`, `results/`, `figures/`) are only readable through an explicit
`CDK2_CAMPAIGN=historical`. Inputs are frozen copies of the cached ChEMBL download and the 4KD1 files, with hashes equal to the historical ones; each stage writes a manifest with the
hashes of the files it consumed. A driver (`scripts/campaign_run.py`) runs the 69 stages from a dependency table and keeps a ledger (planned / running / completed / failed).
Historical files were hashed before the first campaign stage (`historical_record/`) and verified again after the last stage: 132 files, 0 changed, 0 missing (after the documents were rewritten, only README.md, CLAUDE.md and the append-only decisions.md differ; the first two are preserved in `historical_record/`).

**Why.** The earlier results could not be told apart from the output of repaired code, and several defects (docking job/pose association, order-dependent counts, candidates scored without
the standardization applied to training data) made some earlier numbers hard to interpret. A new directory lets old and new coexist without either being edited.

**Alternatives considered.** Overwrite in place (rejected: destroys the evidence); version suffixes on file names (rejected: easy to read the wrong file).

**Cost.** Wall-clock times are in the ledger (docking alone took about 5 h per policy scope for ~2,000 jobs); an extra layer of machinery that a reader must learn (see the README).

**How I'd defend it.** "I wanted the old numbers to stay exactly as they were and the new numbers to be impossible to confuse with them; the price is a driver and some environment variables."

**What would change my mind.** A simpler way to get the same guarantee, e.g. one Git tag per run with read-only results.

## D-41: Frozen plan; the two preparation policies are evaluated as one combined intervention; H3 has a validity gate

**Date:** 2026-10-03 (plan frozen before any stage ran); amendments 2026-10-04
**Stage:** 5, 8, 10

**Decision.** The execution and analysis plan (`campaigns/c2026-10-03/plan/CAMPAIGN_PLAN.md`, hash recorded in the ledger) fixed, before any output existed: the unit of replication (GA run,
n = 5; split for the five scaffold-split surrogates), distinct-molecule metrics as primary, cluster bootstrap for docking statistics, operational definitions of H1-H5, and a **validity gate for H3**:
a held-out comparator can only show "no transfer" or "transfer" through Vina if Vina first separates held-out actives from property-matched decoys (lower bound of the 95% cluster-bootstrap
interval of AUC above 0.5, raw and stratified by heavy-atom count). If it fails, H3 is "not assessable" and the docking-based H5 contingency is not run. The legacy / corrected
policy change bundles three things (standardize like training data, closed-shell only, deduplicate after standardization) and is evaluated only as that bundle: no component ablation.

**Why.** The historical docking result ("no discrimination, H3 inconclusive") was seen before the plan was written; a gate decided in advance stops "Vina did not discriminate" from being
reinterpreted as "optimized molecules did not transfer" after the fact. Pre-specifying n and the unit stops pseudo-replication (generated molecules from one GA run are not independent).

**Alternatives considered.** Report P_sup against held-out actives regardless of the gate (rejected: uninterpretable if Vina carries no activity signal); ablate the three components (rejected: out of budget and not asked).

**Cost.** H3, the most decision-relevant hypothesis, ends as "not assessable" (see D-43). With n = 5 the smallest possible two-sided exact sign-flip p is 0.0625, so p-values are secondary.

**How I'd defend it.** "I decided what would count as evidence before I had any, and I wrote down that Vina had to prove it could tell actives from decoys first. It could not."

**What would change my mind.** A different orthogonal signal (a different docking program, a ligand-based model not trained on this data, or measured data) that passes a gate of the same form.

## D-42: Repairs made while the campaign ran (none touched a seed, parameter or threshold)

**Date:** 2026-10-04
**Stage:** 8, 11, 12

**Decision.** Five defects were fixed and the affected stages rerun in full: `05g` read renamed metric keys and unpacked two values from a three-value function (`holdout_eval`, both policies);
the pose validator compared a numeric-looking `job_id` after RDKit converted it to a float, refusing two correct poses (`legacy/dock` export); the `figures` stage declared file names the script does not
write; figure titles ran off the page; `11b_table1.py` still expected the pre-repair file names and sets. See `plan/CAMPAIGN_PLAN_AMENDMENTS.md` (A4-A7). Failed attempts stay in the ledger and `superseded/`.

**Why.** All five are implementation defects visible as crashes or defective artifacts, not as unwelcome results. The pose-validator case is the instructive one: a stricter integrity check
produced a false alarm, and the response was to prove the poses correct (identical InChIKeys) and fix the check, not to loosen it.

**Alternatives considered.** Bypass the validator for the two jobs (rejected: defeats the check); delete and redock them (rejected: nothing was wrong with them).

**Cost.** About 2.5 h of idle docking time overnight (the driver stopped at the failed export); the code identity now differs between stages (stage manifests record it).

**How I'd defend it.** "Every fix was to code that crashed or drew something wrong; the experiment's seeds, parameters and rules were frozen and never touched, and the failed runs are on record."

**What would change my mind.** Evidence that any fix changed what a stage computes (rerun and compare; the regression tests cover the validator and the metric keys).

## D-43: Campaign outcome: H3 is not assessable; the docking-based H5 contingency was not run

**Date:** 2026-10-05
**Stage:** 8, 9, 10

**Decision.** Under the pre-registered gate, Vina does not discriminate held-out actives from matched decoys in either policy: AUC (actives vs decoys, -Vina) with 95% cluster-bootstrap interval, raw /
size-stratified, legacy M 0.59 [0.46, 0.72] / 0.61 [0.48, 0.76], S 0.57 [0.43, 0.77] / 0.57 [0.42, 0.77]; corrected M 0.55 [0.43, 0.69] / 0.57 [0.46, 0.72], S 0.54 [0.39, 0.77] / 0.53 [0.36, 0.78]. All lower bounds
are below 0.5, so H3 is "not assessable with Vina on this chemistry" for both groups and both policies, the pre-registered H5 docking contingency is not triggered, and P_sup values are tabulated
as descriptive only. The in-sample group P fails the gate as well (0.57 [0.46, 0.68] legacy; 0.60 [0.49, 0.70] corrected). Spearman between -Vina and measured pActivity over reference sets is 0.01-0.10 with intervals including zero.

**Why.** This is what the numbers say under the rule written before they existed. The effect of heavy-atom count on Vina (-0.04 to -0.07 kcal/mol per heavy atom, Spearman -0.3 to -0.5) is larger than any activity signal measured here.

**Alternatives considered.** Run the H5 docking contingency anyway (rejected: breaks the pre-registered rule and would produce numbers without a meaning); switch to a different docking program (non-goal and a new experiment).

**Cost.** H3 and the transfer half of H5 cannot be answered by this project. They are reported as unanswered, not as negative.

**How I'd defend it.** "A docking score that cannot tell known actives from decoys cannot tell me whether generated molecules transfer; saying so is the result."

**What would change my mind.** A gate pass with a different or better-validated scoring signal.

## D-44: Interpretation of H1-H5 from the campaign (policy-by-policy; the legacy arm is kept)

**Date:** 2026-10-05
**Stage:** 5, 10, 11

**Decision.** (i) H1 holds as "rises" for multi_real and activity_only in both policies (net rise of the median real prediction +1.07 / +1.91 legacy, +1.13 / +1.93 corrected; 5/5 runs; Spearman >= 0.9 in 5/5), but the literal "monotonically" holds only for
activity_only (5/5 runs with every step non-decreasing; multi_real 0/5, 86% of steps non-decreasing). (ii) H2 as worded (similarity falls as prediction climbs) holds for multi_real in 4 of 5 runs in each policy (median Spearman -0.95 / -0.92),
but it is not attributable to the activity objective: druglike_only, which has no activity term, ends at the same similarity (0.36), and activity_only moves the other way (final median similarity ~0.70; the literal criterion fails: 0 of 5 runs with a negative Spearman). (iii) H4 holds as
worded: the scrambled-surrogate arm shows Spearman(generation, followed prediction) >= 0.9 in 5/5 runs, with a smaller standardized rise (ratio ~0.7 to the real arm), and its molecules are no better than drug-likeness-only molecules under the real
surrogate (median difference -0.02 legacy, -0.13 corrected) while multi_real's are +1.6 higher. (iv) H5: attenuation of predicted activity (negative paired difference in at least 4 of 5 seeds) appears only for activity_only and only at some floors (legacy: molecule 0.4/0.6/0.7, primary 0.6, scaffold 0.5/0.6/0.7; corrected: molecule 0.4/0.7, primary 0.7, scaffold 0.6/0.7), never for multi_real, whose constrained runs
end at about equal or higher median predicted activity than their unconstrained twins (-0.01 to +0.46); the exchange rate has an interval excluding 0 for scaffold-holdout activity_only at floor 0.7 in both policies (and at 0.4-0.6 corrected), and for molecule-holdout only at 0.7 corrected with a negative sign; elsewhere the point estimate is undefined or the interval spans 0 (a few resample intervals are printed where the point estimate is undefined; read `13_h5_summary.csv`, column `exchange_defined_share_of_resamples`, before using any of them). The legacy arm is retained next to the corrected arm in every table; conclusions of (i)-(iv) are the same in both policies.

**Why.** These are the endpoints as defined in the frozen plan. Where the wording of a hypothesis is stronger than the data (monotonic) or answers a different question than the one that matters (similarity falls vs similarity falls because of the activity objective) both versions are reported.

**Alternatives considered.** Report only the corrected policy (rejected: hides the comparison the campaign exists to make); smooth the "monotonic" claim into "rises" (rejected: the plan separated them).

**Cost.** A longer, less headline-friendly account.

**How I'd defend it.** "I report each hypothesis as written and, where it passes only in a weaker form, say which form and why the stronger one fails."

**What would change my mind.** A multi-start or longer GA where multi_real's steps are monotone; or a similarity measure that is size-independent and still shows attribution.

## D-45: Documentation-and-archive cleanup after the independent review (no science changed)

**Date:** 2026-10-05
**Stage:** 11 (documentation), repository organisation

**Decision.** After the independent review (`review-2026-10-05/REVIEW.md`, `FILE_RETENTION.csv`) the project was reorganised in a pass limited to documentation and file moves: nothing was run (no pipeline stage, test, analysis, figure or experiment), no scientific code or number was changed, and the frozen campaign (`campaigns/c2026-10-03/`) and `historical_record/` were not touched.
(1) `README.md` rewritten as a concise entry point; `CLAUDE.md` made consistent with it (original hypotheses kept as history; the old stage table labelled historical); `KNOWN_ISSUES.md` added; pre-edit copies kept in `archive/documentation_before_cleanup/`.
(2) 123 historical files (top-level `results/`, `figures/`, `data/` except the three campaign inputs) moved intact to `archive/historical_pre_campaign/` with original relative paths; SHA-256 recorded before and after for each (`archive/archive_manifest.csv`, all equal); restoration without overwriting described in `archive/RESTORE.md` (script written, not run).
(3) `presentation/` created with copies of selected existing figures, captions that state quantity, policy, counting and limitations, and a claim-to-evidence table.
Statements corrected from the earlier README/HANDOFF text, against the saved tables and code: the H5 exchange-rate sentence (the saved table has intervals of both signs and conflicts with it; the estimator is unreliable, so exchange rates are provisional and carry no headline claim); endpoint statistics are population-weighted raw-string statistics, not the standardized-distinct view the plan names as primary, and docking sets deduplicate raw SMILES;
"transfer" is now always named (surrogate prediction, fingerprint proximity, docking score); an inconclusive docking benchmark is not described as proof that Vina has no signal; "docking avoids clashes by construction" and any suggestion that a PoseBusters pass shows binding were removed; "standard library only" was corrected (unittest is the runner, the tests need the scientific packages); the plan is described as frozen before this campaign but written after the earlier exploration;
saved test results are attributed to the review's run of 2026-10-05, not to this cleanup. `campaigns/c2026-10-03/HANDOFF.md` still contains superseded wording; it was not edited because it is inside the frozen campaign (KNOWN_ISSUES #10).

**Why.** The review found statements in the README that the saved tables contradict and a counting definition the code does not implement. Fixing the wording does not require, and must not be mistaken for, repairing the analyses: those need new computation and would be a post-campaign amendment.

**Alternatives considered.** Repair the endpoint estimator and recompute (rejected for this pass by instruction, and it would change results); delete the historical files or caches to tidy the tree (rejected: evidence; the three campaign input paths and everything the campaign reads stay in place); edit HANDOFF.md (rejected: inside the frozen, inventoried campaign).

**Cost.** `verify-historical` now reports the 123 archived files as missing until restored; the tests that read historical files skip themselves; HANDOFF.md and README now differ.

**How I'd defend it.** "I corrected the claims that the saved tables did not support, labelled what I could not verify as provisional, and moved the old results out of the way without losing a byte; I did not pretend that rewording fixes the estimator."

**What would change my mind.** A recomputation of the H5 and counting analyses with a corrected estimator and standardized-distinct views: the documents would then be rewritten from the new tables as a recorded amendment.

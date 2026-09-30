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

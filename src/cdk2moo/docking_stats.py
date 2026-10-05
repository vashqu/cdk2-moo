"""
Statistics for the docking analysis, as defined in the campaign plan (section 3).

Conventions. A Vina score is "better" when MORE NEGATIVE. P_sup(a, b) is the probability that a random molecule of sample a scores better than a
random molecule of sample b (ties count half); 0.5 means no difference. AUC(actives vs reference) uses -Vina as the score, so it equals P_sup(actives, reference).

Uncertainty. Molecules are not independent replicates: generated molecules and decoys descend from a few GA runs, and reference molecules cluster in
chemical series. The cluster bootstrap resamples whole CLUSTERS with replacement, separately in the two samples compared (clusters = GA seed for generated
molecules and decoys, Murcko scaffold for reference molecules), and recomputes the statistic; intervals are 95% percentiles of 2,000 resamples (generator seed 42).

Size. Vina scores improve with heavy-atom count, so statistics are also computed within heavy-atom strata and combined with weights equal to the number of
molecule pairs in each stratum. A stratum needs at least MIN_PER_STRATUM molecules in each sample to contribute.
"""

import numpy as np

STRATA = [(20, 24), (25, 29), (30, 34), (35, 39)]
MIN_PER_STRATUM = 5
N_BOOT = 2000


def p_better(scores_a, scores_b):
    """P(a random molecule of a scores better (lower) than a random molecule of b); ties count half. NaN if a sample is empty."""
    a, b = np.asarray(scores_a, float), np.asarray(scores_b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    less = (a[:, None] < b[None, :]).mean()
    ties = (a[:, None] == b[None, :]).mean()
    return float(less + 0.5 * ties)


def stratum_of(heavy_atoms):
    """Label (index into STRATA) of each molecule, -1 outside every stratum."""
    atoms = np.asarray(heavy_atoms)
    label = np.full(len(atoms), -1)
    for k, (low, high) in enumerate(STRATA):
        label[(atoms >= low) & (atoms <= high)] = k
    return label


def stratified_p_better(scores_a, atoms_a, scores_b, atoms_b):
    """
    P_sup computed inside each heavy-atom stratum (needing MIN_PER_STRATUM molecules in both samples) and averaged with weights n_a * n_b.
    Returns (value, pairs_used, pairs_total); value is NaN if no stratum qualifies.
    """
    a, b = np.asarray(scores_a, float), np.asarray(scores_b, float)
    sa, sb = stratum_of(atoms_a), stratum_of(atoms_b)
    total = len(a) * len(b)
    num, den = 0.0, 0.0
    for k in range(len(STRATA)):
        ia, ib = a[sa == k], b[sb == k]
        if len(ia) >= MIN_PER_STRATUM and len(ib) >= MIN_PER_STRATUM:
            weight = len(ia) * len(ib)
            num += weight * p_better(ia, ib)
            den += weight
    return (num / den if den else float("nan")), int(den), int(total)


def _group_by_cluster(clusters):
    clusters = np.asarray(clusters, dtype=object)
    labels, inverse = np.unique(clusters.astype(str), return_inverse=True)
    return [np.where(inverse == i)[0] for i in range(len(labels))]


def cluster_bootstrap(statistic, a, clusters_a, b, clusters_b, n_boot=N_BOOT, seed=42):
    """
    95% cluster-bootstrap interval of statistic(a_rows, b_rows). `a` and `b` are arrays of row indices' payloads: this function resamples row INDICES of
    each sample by whole clusters and calls statistic(index_a, index_b). Returns (low, high, n_valid) over resamples with a defined value.
    """
    rng = np.random.default_rng(seed)
    groups_a, groups_b = _group_by_cluster(clusters_a), _group_by_cluster(clusters_b)
    values = []
    for _ in range(n_boot):
        ia = np.concatenate([groups_a[i] for i in rng.integers(0, len(groups_a), len(groups_a))]) if groups_a else np.array([], int)
        ib = np.concatenate([groups_b[i] for i in rng.integers(0, len(groups_b), len(groups_b))]) if groups_b else np.array([], int)
        value = statistic(ia, ib)
        if np.isfinite(value):
            values.append(value)
    if not values:
        return float("nan"), float("nan"), 0
    return float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5)), len(values)


def compare_samples(score_a, atoms_a, clusters_a, score_b, atoms_b, clusters_b, n_boot=N_BOOT, seed=42):
    """
    P_sup of sample a over sample b, raw and size-stratified, each with a cluster-bootstrap interval. Returns a dict.
    """
    score_a, score_b = np.asarray(score_a, float), np.asarray(score_b, float)
    atoms_a, atoms_b = np.asarray(atoms_a), np.asarray(atoms_b)
    out = {"n_a": len(score_a), "n_b": len(score_b), "clusters_a": len(set(map(str, clusters_a))), "clusters_b": len(set(map(str, clusters_b)))}
    out["p_sup"] = p_better(score_a, score_b)
    out["p_sup_stratified"], out["pairs_used"], out["pairs_total"] = stratified_p_better(score_a, atoms_a, score_b, atoms_b)
    if len(score_a) and len(score_b):
        raw = lambda ia, ib: p_better(score_a[ia], score_b[ib])                                       # noqa: E731
        strat = lambda ia, ib: stratified_p_better(score_a[ia], atoms_a[ia], score_b[ib], atoms_b[ib])[0]   # noqa: E731
        out["p_sup_low"], out["p_sup_high"], _ = cluster_bootstrap(raw, score_a, clusters_a, score_b, clusters_b, n_boot, seed)
        out["p_sup_strat_low"], out["p_sup_strat_high"], _ = cluster_bootstrap(strat, score_a, clusters_a, score_b, clusters_b, n_boot, seed)
    else:
        out.update({k: float("nan") for k in ("p_sup_low", "p_sup_high", "p_sup_strat_low", "p_sup_strat_high")})
    return out


def validity_gate(auc_raw_low, auc_strat_low, auc_raw, auc_strat):
    """
    The plan's gate for using Vina as an orthogonal signal: the lower bounds of the raw and the stratified AUC (actives vs matched decoys) must exceed
    0.5. Returns (passed, label): label "weak" if the point AUC is below 0.65 in either version, "moderate" otherwise; "failed" or "undefined" if not passed.
    """
    values = [auc_raw_low, auc_strat_low, auc_raw, auc_strat]
    if not all(np.isfinite(values)):
        return False, "undefined"
    if auc_raw_low > 0.5 and auc_strat_low > 0.5:
        return True, "moderate" if min(auc_raw, auc_strat) >= 0.65 else "weak"
    return False, "failed"

"""
Per-run trajectory statistics for the hypothesis endpoints defined in the campaign plan (section 3).

A RUN is one arm x one GA seed. For a run, the population at each generation is a table of molecules. Statistics about generated chemistry use only
molecules born after generation 0; generation 0 (the starting training molecules) is shown as the reference. A generation counts as QUALIFIED for step
and trend statistics if it has at least MIN_GENERATED generated molecules (generation 0 always qualifies); unqualified generations are reported, not hidden.
"""

import numpy as np
from scipy.stats import spearmanr

MIN_GENERATED = 10


def generation_table(run, column):
    """
    Median of `column` per generation for one run (a DataFrame with generation, birth_generation and `column`).
    Columns returned: generation, n_generated (molecules counted), median, qualified.
    """
    rows = []
    for generation, group in run.groupby("generation"):
        used = group if generation == 0 else group[group["birth_generation"] > 0]
        rows.append({"generation": int(generation), "n_generated": int(len(used)),
                     "median": float(used[column].median()) if len(used) else float("nan"),
                     "qualified": bool(generation == 0 or len(used) >= MIN_GENERATED)})
    return rows


def _qualified(rows, first_generation=0):
    keep = [r for r in rows if r["qualified"] and r["generation"] >= first_generation and np.isfinite(r["median"])]
    return np.array([r["generation"] for r in keep]), np.array([r["median"] for r in keep])


def trend_statistics(rows):
    """
    For a per-generation median series: net rise (last qualified minus generation 0), the fraction of consecutive qualified generation pairs
    with a non-decreasing median (the literal "monotonic" test), the Spearman correlation of generation with the median, and how many generations qualified.
    """
    gens, values = _qualified(rows)
    if len(gens) < 3:
        return {"net_rise": float("nan"), "fraction_nondecreasing": float("nan"), "spearman_generation": float("nan"), "n_qualified": int(len(gens))}
    steps = np.diff(values)
    return {"net_rise": float(values[-1] - values[0]), "fraction_nondecreasing": float(np.mean(steps >= 0)),
            "spearman_generation": float(spearmanr(gens, values)[0]), "n_qualified": int(len(gens))}


def association_across_generations(rows_a, rows_b, first_generation=1):
    """Spearman correlation between two per-generation median series over the generations qualified in BOTH (from `first_generation`)."""
    a = {r["generation"]: r["median"] for r in rows_a if r["qualified"] and r["generation"] >= first_generation and np.isfinite(r["median"])}
    b = {r["generation"]: r["median"] for r in rows_b if r["qualified"] and r["generation"] >= first_generation and np.isfinite(r["median"])}
    common = sorted(set(a) & set(b))
    if len(common) < 3:
        return float("nan"), len(common)
    return float(spearmanr([a[g] for g in common], [b[g] for g in common])[0]), len(common)


def bootstrap_ratio(numerator, denominator, n_boot=2000, seed=42):
    """
    Bootstrap over runs (seeds) of mean(numerator) / (-mean(denominator)): the exchange rate, defined only when the mean denominator is negative
    (something was given up). Returns (estimate, low, high, share_of_resamples_defined); estimate NaN if the point estimate is undefined.
    """
    num, den = np.asarray(numerator, float), np.asarray(denominator, float)
    point = np.nan if den.mean() >= 0 else float(num.mean() / -den.mean())
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(num), len(num))
        if den[idx].mean() < 0:
            values.append(num[idx].mean() / -den[idx].mean())
    if not values:
        return point, float("nan"), float("nan"), 0.0
    return point, float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5)), len(values) / n_boot

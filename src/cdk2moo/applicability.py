"""
Applicability-domain audit: how does the surrogate's error depend on how close a
molecule is to the training set?

Two distance measures, kept separate and compared:
  raw similarity   highest ECFP4 Tanimoto to the training set
  size percentile  that value as a percentile among size-matched training
                   molecules (leave-one-out; see similarity.py)

The error estimates come from Stage 4 (out-of-sample predictions for molecules
the forest never saw). A cell of the validation data with too few distinct
molecules is flagged, and no error is claimed for it.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from cdk2moo import config


def add_region(df, column, edges, names, out):
    """Add a categorical column `out` labelling which region `column` falls in."""
    df[out] = pd.cut(df[column], edges, right=False, labels=names)
    return df


def error_table(df, by):
    """
    Error of the forest within each region (or combination of regions).

    `by` is a list of region columns. Reports rows, DISTINCT molecules (the same
    molecule appears in several seeds' test sets, so rows overstate independent
    evidence), mean absolute error, RMSE and bias (prediction minus measured;
    positive = over-prediction). `supported` is False below config.AD_MIN_UNIQUE
    distinct molecules.
    """
    rows = []
    for key, g in df.groupby(by, observed=True):
        err = g["pred_rf"] - g["y_true"]
        key = key if isinstance(key, tuple) else (key,)
        rows.append({**dict(zip(by, key)),
                     "n_rows": len(g), "n_unique": g["mol"].nunique(),
                     "mae": err.abs().mean(), "rmse": float(np.sqrt((err ** 2).mean())),
                     "bias": err.mean(),
                     "supported": g["mol"].nunique() >= config.AD_MIN_UNIQUE})
    return pd.DataFrame(rows)


def large_error_auc(df, score_col):
    """
    How well does a similarity measure flag molecules the forest gets badly wrong?

    AUC = probability that a randomly chosen large-error molecule has LOWER
    similarity than a randomly chosen small-error one. 0.5 = no information,
    1.0 = perfect. Large error means |error| >= config.AD_LARGE_ERROR.
    """
    bad = (df["pred_rf"] - df["y_true"]).abs() >= config.AD_LARGE_ERROR
    if bad.all() or not bad.any():
        return float("nan")
    return float(roc_auc_score(bad, -df[score_col]))

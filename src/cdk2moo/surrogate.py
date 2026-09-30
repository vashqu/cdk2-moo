"""
The activity surrogate: a random forest that maps a fingerprint to a predicted
pActivity, plus the controls that say what its numbers mean.

The surrogate is an INSTRUMENT, not a result. Later the genetic algorithm will
treat its output as the thing to maximise, and the project asks what happens
when it does. So we do not tune it to look good; we fix its settings in config
and measure how far it can be trusted, with a control for every metric.

Random forest in one paragraph: train many decision trees, each on a random
resample of the molecules and considering a random subset of fingerprint bits
at each split; the prediction is the average of the trees. A tree can only
predict values it has seen in training, so a forest cannot predict outside the
range of its training labels, and for a molecule unlike anything in training the
trees fall back on whatever bits happen to match. That is where extrapolation
error will come from.
"""

import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

from cdk2moo import config


def fit_forest(X, y, seed, max_features=None):
    """Fit the project's random forest. Settings come from config, never tuned."""
    if max_features is None:
        max_features = config.RF_MAX_FEATURES
    forest = RandomForestRegressor(
        n_estimators=config.RF_N_TREES,
        max_features=max_features,
        random_state=seed,
        n_jobs=-1,
    )
    return forest.fit(X, y)


def scramble_labels(y, seed):
    """
    Shuffle the training labels across the training molecules.

    The scrambled labels keep exactly the same distribution of values but break
    every link between structure and activity. A model trained on them has
    nothing real to learn. If it still scores well on the test set, the
    evaluation is broken; if it scores ~0, as it should, it becomes the control
    arm for the GA (hypothesis H4): an optimizer following a surrogate that
    knows nothing, to see what an "optimization curve" looks like by itself.
    """
    return np.random.default_rng(seed).permutation(y)


def regression_metrics(y_true, y_pred):
    """
    R2, RMSE, MAE and Spearman rank correlation of predictions vs measurements.

    R2 is 1 for perfect predictions, 0 for "always predict the TEST-set mean",
    and NEGATIVE if worse than that. (The test-set mean is not knowable in
    advance: when test labels are shifted relative to train, even predicting the
    training mean scores slightly below 0. That is why the "mean" control exists.) Spearman only asks whether the
    ranking is right, which is what an optimizer needs from a surrogate, and is
    insensitive to a constant offset. It is NaN when predictions are constant.
    """
    err = y_pred - y_true
    ss_res = float(np.sum(err ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    if np.ptp(y_pred) == 0:
        rho = float("nan")
    else:
        rho = float(spearmanr(y_true, y_pred)[0])
    return {
        "r2": 1 - ss_res / ss_tot,
        "rmse": float(np.sqrt(np.mean(err ** 2))),
        "mae": float(np.mean(np.abs(err))),
        "spearman": rho,
    }

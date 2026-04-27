"""
Simulate a univariate series and score a deterministic roerich ChangePointDetectionClassifier-style detector.

This uses the roerich window/classification pipeline with a deterministic QDA backend
and a shared stratified split generator so Python and Fortran can be compared exactly.
"""

from __future__ import annotations

import time

import numpy as np
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.preprocessing import StandardScaler


def simulate_roerich_classifier_signal(
    n: int,
    true_cp: int,
    means: tuple[float, float],
    sds: tuple[float, float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with one mean change.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    signal[:true_cp] = rng.normal(means[0], sds[0], size=true_cp)
    signal[true_cp:] = rng.normal(means[1], sds[1], size=n - true_cp)
    return signal


def autoregression_matrix_1d(signal: np.ndarray, periods: int) -> np.ndarray:
    """
    Build the univariate autoregression feature matrix used by roerich.
    """
    n = signal.shape[0]
    x_auto = np.zeros((n, periods), dtype=float)
    for lag in range(periods):
        x_auto[lag:, lag] = signal[: n - lag]
    return x_auto


def park_miller_next(state: int) -> int:
    """
    Advance the Park-Miller generator by one step.
    """
    a = 16807
    m = 2147483647
    q = 127773
    r = 2836
    hi = state // q
    lo = state % q
    test = a * lo - r * hi
    return test if test > 0 else test + m


def lcg_permutation(n: int, seed: int) -> np.ndarray:
    """
    Generate the same Fisher-Yates permutation as the Fortran helper.
    """
    state = abs(seed) % 2147483646 + 1
    perm = np.arange(n, dtype=int)
    for i in range(n - 1, 0, -1):
        state = park_miller_next(state)
        j = (state - 1) % (i + 1)
        perm[i], perm[j] = perm[j], perm[i]
    return perm


def klsym_metric(ref_preds: np.ndarray, test_preds: np.ndarray) -> float:
    """
    Symmetric KL score used by roerich's classifier detector.
    """
    return float(
        np.mean(np.log(test_preds + 1e-3))
        - np.mean(np.log(1.0 - test_preds + 1e-3))
        + np.mean(np.log(1.0 - ref_preds + 1e-3))
        - np.mean(np.log(ref_preds + 1e-3))
    )


def cpdc_qda_score(x_ref: np.ndarray, x_test: np.ndarray, seed: int) -> float:
    """
    Deterministic QDA/KL_sym score for one reference/test pair.
    """
    x = np.vstack((x_ref, x_test))
    y = np.hstack((np.zeros(len(x_ref), dtype=int), np.ones(len(x_test), dtype=int)))
    x = StandardScaler().fit_transform(x)

    idx0 = np.flatnonzero(y == 0)
    idx1 = np.flatnonzero(y == 1)
    idx0 = idx0[lcg_permutation(idx0.size, seed + 37)]
    idx1 = idx1[lcg_permutation(idx1.size, seed + 74)]

    ntest0 = idx0.size // 2
    ntest1 = idx1.size // 2
    test_idx = np.concatenate((idx0[:ntest0], idx1[:ntest1]))
    train_idx = np.concatenate((idx0[ntest0:], idx1[ntest1:]))

    clf = QuadraticDiscriminantAnalysis(store_covariance=True, reg_param=0.01)
    clf.fit(x[train_idx], y[train_idx])
    y_pred = clf.predict_proba(x[test_idx])[:, 1]
    return klsym_metric(y_pred[y[test_idx] == 0], y_pred[y[test_idx] == 1])


def fit_roerich_classifier(
    signal: np.ndarray,
    window_size: int,
    periods: int = 1,
    step: int = 1,
    seed: int = 2357,
) -> tuple[int, float]:
    """
    Run the deterministic raw ChangePointDetectionClassifier-style score.
    """
    x_auto = autoregression_matrix_1d(signal, periods)
    best_cp = -1
    best_score = -np.inf
    for t_idx in range(2 * window_size, signal.shape[0] + 1, step):
        x_ref = x_auto[t_idx - 2 * window_size : t_idx - window_size]
        x_test = x_auto[t_idx - window_size : t_idx]
        score = cpdc_qda_score(x_ref, x_test, seed + t_idx)
        if score > best_score:
            best_score = score
            best_cp = (t_idx - 1) - window_size
    return int(best_cp), float(best_score)


def main() -> None:
    t0 = time.perf_counter()
    n = 400
    true_cp = 200
    window_size = 40
    periods = 1
    step = 1
    signal = simulate_roerich_classifier_signal(n, true_cp, (0.0, 3.0), (1.0, 1.0), 101)
    cp, score = fit_roerich_classifier(signal, window_size, periods, step)

    print("n =", n)
    print("true changepoint =", true_cp)
    print("window_size =", window_size)
    print("periods =", periods)
    print("step =", step)
    print("base_classifier = qda")
    print("metric = klsym")
    print("estimated changepoint =", cp)
    print(f"max raw score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

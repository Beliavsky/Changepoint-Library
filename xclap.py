"""
Simulate a univariate labeled series and score a CLaP-style state classifier.

This mirrors the window-dataset construction from claspy's CLaP workflow, but
uses a deterministic nearest-centroid classifier because the aeon classifier
stack required by claspy.clap is not available in this environment.
"""

from __future__ import annotations

import time
from collections import Counter

import numpy as np


def simulate_clap_signal(
    n: int,
    regime_starts: list[int],
    phis: list[float],
    sigmas: list[float],
    regime_labels: list[int],
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Simulate a univariate AR(1) series with repeated labeled states.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    state_labels = np.empty(n, dtype=int)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        prev = 0.0
        for t in range(start, end):
            prev = phis[i] * prev + rng.normal(loc=0.0, scale=sigmas[i])
            signal[t] = prev
        state_labels[start:end] = regime_labels[i]
    return signal, state_labels


def lcg_permutation(n: int, seed: int) -> np.ndarray:
    """
    Generate the same Fisher-Yates permutation as the Fortran helper.
    """
    a = 16807
    m = 2147483647
    q = 127773
    r = 2836
    state = abs(seed) % (m - 1) + 1
    perm = np.arange(n, dtype=int)
    for i in range(n - 1, 0, -1):
        hi = state // q
        lo = state % q
        test = a * lo - r * hi
        state = test if test > 0 else test + m
        j = (state - 1) % (i + 1)
        perm[i], perm[j] = perm[j], perm[i]
    return perm


def create_clap_dataset(
    signal: np.ndarray,
    state_labels: np.ndarray,
    window_size: int,
    stride: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Create the labeled subsequence dataset used by the CLaP-style classifier.
    """
    cps = np.flatnonzero(state_labels[:-1] != state_labels[1:]) + 1
    starts = np.arange(0, signal.shape[0] - window_size + 1, stride, dtype=int)
    keep = np.ones(starts.shape[0], dtype=bool)

    for cp in cps:
        left = cp - window_size // 2
        right = cp - 1
        keep &= ~((starts >= left) & (starts <= right))

    kept_starts = starts[keep]
    xmat = np.column_stack([signal[start : start + window_size] for start in kept_starts])
    y = state_labels[kept_starts].astype(int, copy=False)
    return xmat, y


def subselect_clap_dataset(
    xmat: np.ndarray,
    y: np.ndarray,
    sample_cap: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Cap the number of windows per label and apply a deterministic shuffle.
    """
    labels = np.unique(y)
    keep: list[int] = []
    for i, label in enumerate(labels, start=1):
        idx = np.flatnonzero(y == label)
        perm = lcg_permutation(idx.shape[0], seed + 37 * i)
        idx = idx[perm]
        take = min(sample_cap, idx.shape[0])
        keep.extend(idx[:take].tolist())
    keep_arr = np.array(keep, dtype=int)
    perm = lcg_permutation(keep_arr.shape[0], seed + 911)
    keep_arr = keep_arr[perm]
    return xmat[:, keep_arr], y[keep_arr]


def crossval_centroid(
    xmat: np.ndarray,
    y: np.ndarray,
    n_splits: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Run deterministic K-fold nearest-centroid classification.
    """
    n_samples = y.shape[0]
    perm = lcg_permutation(n_samples, seed)
    y_true = np.empty(n_samples, dtype=int)
    y_pred = np.empty(n_samples, dtype=int)

    for fold in range(n_splits):
        fold_start = (fold * n_samples) // n_splits
        fold_end = ((fold + 1) * n_samples) // n_splits
        test_mask = np.zeros(n_samples, dtype=bool)
        test_mask[fold_start:fold_end] = True

        test_idx = perm[test_mask]
        train_idx = perm[~test_mask]
        xtrain = xmat[:, train_idx]
        ytrain = y[train_idx]
        labels = np.unique(ytrain)
        centroids = np.column_stack([xtrain[:, ytrain == label].mean(axis=1) for label in labels])

        y_true[test_idx] = y[test_idx]
        for idx in test_idx:
            dists = np.sum((centroids - xmat[:, [idx]]) ** 2, axis=0)
            y_pred[idx] = int(labels[int(np.argmin(dists))])

    return y_true, y_pred


def macro_f1_labels(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Compute macro F1 across the labels present in y_true.
    """
    labels = np.unique(y_true)
    total = 0.0
    for label in labels:
        tp = np.sum((y_true == label) & (y_pred == label))
        fp = np.sum((y_true != label) & (y_pred == label))
        fn = np.sum((y_true == label) & (y_pred != label))
        denom = 2 * tp + fp + fn
        if denom > 0:
            total += 2.0 * tp / denom
    return total / labels.shape[0]


def fit_clap_centroid(
    signal: np.ndarray,
    state_labels: np.ndarray,
    window_size: int,
    n_splits: int,
    sample_cap: int = 1000,
    seed: int = 2357,
) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Score a deterministic nearest-centroid analog of claspy's CLaP pipeline.
    """
    stride = max(1, window_size // 2)
    xmat, y = create_clap_dataset(signal, state_labels, window_size, stride)
    xmat, y = subselect_clap_dataset(xmat, y, sample_cap, seed)
    y_true, y_pred = crossval_centroid(xmat, y, n_splits, seed)
    score = macro_f1_labels(y_true, y_pred)
    return y_true, y_pred, score


def main() -> None:
    t0 = time.perf_counter()
    seed = 211
    n = 480
    window_size = 20
    n_splits = 5
    sample_cap = 1000
    regime_starts = [0, 80, 160, 240, 320, 400]
    phis = [0.8, -0.6, 0.2, 0.8, -0.6, 0.2]
    sigmas = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
    regime_labels = [1, 2, 3, 1, 2, 3]

    signal, state_labels = simulate_clap_signal(
        n, regime_starts, phis, sigmas, regime_labels, seed
    )
    y_true, y_pred, score = fit_clap_centroid(
        signal, state_labels, window_size, n_splits, sample_cap
    )

    print("n =", n)
    print("window_size =", window_size)
    print("n_splits =", n_splits)
    print("sample_cap =", sample_cap)
    print("regime starts =", regime_starts)
    print("regime labels =", regime_labels)
    print("phis =", phis)
    print("sigmas =", sigmas)
    print("true label counts =", dict(sorted(Counter(y_true.tolist()).items())))
    print("pred label counts =", dict(sorted(Counter(y_pred.tolist()).items())))
    print("true labels =", y_true.tolist())
    print("pred labels =", y_pred.tolist())
    print(f"macro_f1 = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

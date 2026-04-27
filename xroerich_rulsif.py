"""
Simulate a univariate series and score a deterministic roerich RuLSIF-style detector.
"""

from __future__ import annotations

import time

import numpy as np


def simulate_roerich_rulsif_signal(
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


def one_side_rulsif_score(
    x_train: np.ndarray,
    x_eval: np.ndarray,
    y_train: np.ndarray,
    y_eval: np.ndarray,
    alpha: float,
    l2: float,
) -> float:
    """
    One-direction linear RuLSIF PE score.
    """
    phi_train = np.column_stack((np.ones(x_train.shape[0]), x_train))
    phi_eval = np.column_stack((np.ones(x_eval.shape[0]), x_eval))
    phi0 = phi_train[y_train == 0]
    phi1 = phi_train[y_train == 1]
    hmat = (1.0 - alpha) / phi0.shape[0] * (phi0.T @ phi0) + alpha / phi1.shape[0] * (phi1.T @ phi1)
    hmat = hmat + l2 * np.eye(hmat.shape[0])
    h = phi1.mean(axis=0)
    theta = np.linalg.solve(hmat, h)
    ratios = phi_eval @ theta
    ratios = np.maximum(ratios, 0.0)
    return 0.5 * ratios[y_eval == 1].mean() - 0.5


def rulsif_linear_pesym_score(
    x_ref: np.ndarray,
    x_test: np.ndarray,
    seed: int,
    alpha: float = 0.05,
    l2: float = 1e-3,
) -> float:
    """
    Deterministic two-sided linear RuLSIF score for one reference/test pair.
    """
    x = np.vstack((x_ref, x_test))
    y = np.hstack((np.zeros(len(x_ref), dtype=int), np.ones(len(x_test), dtype=int)))
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    sd[sd < 1e-12] = 1.0
    x = (x - mu) / sd

    idx0 = np.flatnonzero(y == 0)
    idx1 = np.flatnonzero(y == 1)
    idx0 = idx0[lcg_permutation(idx0.size, seed + 37)]
    idx1 = idx1[lcg_permutation(idx1.size, seed + 74)]
    ntest0 = idx0.size // 2
    ntest1 = idx1.size // 2
    test_idx = np.concatenate((idx0[:ntest0], idx1[:ntest1]))
    train_idx = np.concatenate((idx0[ntest0:], idx1[ntest1:]))

    score_right = one_side_rulsif_score(
        x[train_idx], x[test_idx], y[train_idx], y[test_idx], alpha, l2
    )
    score_left = one_side_rulsif_score(
        x[test_idx], x[train_idx], 1 - y[test_idx], 1 - y[train_idx], alpha, l2
    )
    return float(score_right + score_left)


def fit_roerich_rulsif(
    signal: np.ndarray,
    window_size: int,
    periods: int = 1,
    step: int = 1,
    seed: int = 2357,
    alpha: float = 0.05,
    l2: float = 1e-3,
) -> tuple[int, float]:
    """
    Run the deterministic raw RuLSIF-style score.
    """
    x_auto = autoregression_matrix_1d(signal, periods)
    best_cp = -1
    best_score = -np.inf
    for t_idx in range(2 * window_size, signal.shape[0] + 1, step):
        x_ref = x_auto[t_idx - 2 * window_size : t_idx - window_size]
        x_test = x_auto[t_idx - window_size : t_idx]
        score = rulsif_linear_pesym_score(x_ref, x_test, seed + t_idx, alpha, l2)
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
    alpha = 0.05
    l2 = 1e-3
    signal = simulate_roerich_rulsif_signal(n, true_cp, (0.0, 3.0), (1.0, 1.0), 101)
    cp, score = fit_roerich_rulsif(signal, window_size, periods, step, 2357, alpha, l2)

    print("n =", n)
    print("true changepoint =", true_cp)
    print("window_size =", window_size)
    print("periods =", periods)
    print("step =", step)
    print("base_regressor = linear_rulsif")
    print("metric = pesym")
    print("alpha =", alpha)
    print("l2 =", l2)
    print("estimated changepoint =", cp)
    print(f"max raw score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

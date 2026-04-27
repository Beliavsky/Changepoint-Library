"""
Compare several ruptures cost models on a univariate series with changes
in both mean and variance.
"""

from __future__ import annotations

import time

import numpy as np
import ruptures as rpt


def simulate_mean_variance_signal(
    n: int,
    bkps: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise mean and variance.
    """
    if len(bkps) != len(means) or len(bkps) != len(sds):
        raise ValueError("bkps, means, and sds must have the same length")
    if bkps[-1] != n:
        raise ValueError("final breakpoint must equal n")

    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    start = 0
    for end, mean, sd in zip(bkps, means, sds):
        signal[start:end] = rng.normal(loc=mean, scale=sd, size=end - start)
        start = end
    return signal


def fit_cost_model(signal: np.ndarray, model: str, n_bkps: int, min_size: int) -> tuple[list[int], float]:
    """
    Fit Dynp with one ruptures cost model and return breakpoints and elapsed time.
    """
    t0 = time.perf_counter()
    algo = rpt.Dynp(model=model, min_size=min_size, jump=1).fit(signal.reshape(-1, 1))
    bkps = algo.predict(n_bkps=n_bkps)
    elapsed = time.perf_counter() - t0
    return bkps, elapsed


def main() -> None:
    t0 = time.perf_counter()

    seed = 11
    n = 600
    true_bkps = [150, 300, 450, n]
    means = [0.0, 2.5, -1.0, 1.5]
    sds = [1.0, 2.2, 0.6, 1.8]
    n_bkps = len(true_bkps) - 1
    min_size = 30

    signal = simulate_mean_variance_signal(n, true_bkps, means, sds, seed)

    models = [
        ("l1", "robust location shifts"),
        ("l2", "mean shifts"),
        ("normal", "mean and variance shifts"),
        ("rbf", "general distribution shifts"),
    ]

    print("n =", n)
    print("true bkps =", true_bkps)
    print("segment means =", means)
    print("segment sds =", sds)
    print("n_bkps =", n_bkps)
    print("min_size =", min_size)
    print()

    for model, note in models:
        bkps, elapsed = fit_cost_model(signal, model, n_bkps, min_size)
        print(f"model = {model}")
        print("note =", note)
        print("estimated bkps =", bkps)
        print(f"elapsed seconds = {elapsed:.3f}")
        print()

    print(f"wall time elapsed (s) = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

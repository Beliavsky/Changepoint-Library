"""
Read a univariate series from a text file and compare several ruptures cost models
for changes in both mean and variance.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import ruptures as rpt


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
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcosts_mean_variance_data.txt"
    n_bkps = 3
    min_size = 30

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n = signal.shape[0]
    models = [
        ("l1", "robust location shifts"),
        ("l2", "mean shifts"),
        ("normal", "mean and variance shifts"),
        ("rbf", "general distribution shifts"),
    ]

    print("file =", data_file)
    print("n =", n)
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

"""
Simulate a univariate series and score it with a deterministic Gaussian-kernel ULSIF detector.
"""

from __future__ import annotations

import time

from xrulsif_utils import fit_rulsif_gaussian, simulate_rulsif_signal


def main() -> None:
    t0 = time.perf_counter()
    n = 400
    true_cp = 200
    window_length = 10
    n_windows = 50
    lag = 50
    sigma = 2.0
    lambda_reg = 1e-2
    scoring_step = 1

    signal = simulate_rulsif_signal(n, true_cp, (0.0, 3.0), (1.0, 1.0), 101)
    _, cp, best_score = fit_rulsif_gaussian(
        signal,
        window_length=window_length,
        n_windows=n_windows,
        lag=lag,
        alpha=0.0,
        sigma=sigma,
        lambda_reg=lambda_reg,
        scoring_step=scoring_step,
        symmetric=True,
    )

    print("n =", n)
    print("true changepoint =", true_cp)
    print("window_length =", window_length)
    print("n_windows =", n_windows)
    print("lag =", lag)
    print("alpha = 0.0")
    print("sigma =", sigma)
    print("lambda =", lambda_reg)
    print("scoring_step =", scoring_step)
    print("estimated changepoint =", cp)
    print(f"max raw score = {best_score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

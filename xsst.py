"""
Simulate a frequency-change signal and score it with changepoynt SST.
"""

from __future__ import annotations

import time

from xsst_utils import fit_sst_naive, simulate_frequency_change_signal


def main() -> None:
    t0 = time.perf_counter()
    n_per_segment = 320
    true_cp = n_per_segment
    period_before = 48
    period_after = 14
    noise = 0.02
    window_length = 48
    n_windows = 48
    lag = 16
    rank = 2
    scoring_step = 1

    signal, _ = simulate_frequency_change_signal(
        n_per_segment=n_per_segment,
        period_before=period_before,
        period_after=period_after,
        noise=noise,
        seed=5678,
    )
    _, cp, best_score = fit_sst_naive(
        signal,
        window_length=window_length,
        n_windows=n_windows,
        lag=lag,
        rank=rank,
        scoring_step=scoring_step,
    )

    print("n =", signal.size)
    print("true changepoint =", true_cp)
    print("window_length =", window_length)
    print("n_windows =", n_windows)
    print("lag =", lag)
    print("rank =", rank)
    print("scoring_step =", scoring_step)
    print("method = naive")
    print("estimated changepoint =", cp)
    print(f"max raw score = {best_score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

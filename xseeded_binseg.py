"""
Simulate a univariate Gaussian series and detect changepoints with skchange.SeededBinarySegmentation.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from skchange.change_detectors import SeededBinarySegmentation


def simulate_seeded_binseg_signal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise mean and variance.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.normal(loc=means[i], scale=sds[i], size=end - start)
    return signal


def fit_seeded_binseg(signal: np.ndarray) -> list[int]:
    """
    Run SeededBinarySegmentation and return detected changepoints.
    """
    n = signal.shape[0]
    penalty = 2.0 * np.log(n)
    detector = SeededBinarySegmentation(
        penalty=penalty,
        max_interval_length=max(200, n // 8),
        growth_factor=1.5,
        selection_method="greedy",
    )
    result = detector.fit_predict(pd.DataFrame({"x": signal}))
    return result["ilocs"].to_list()


def main() -> None:
    t0 = time.perf_counter()
    seed = 61
    n = 20_000
    true_cps = [0, 4_000, 8_000, 12_000, 16_000]
    means = [0.0, 4.0, -2.0, 3.0, -1.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]
    penalty = 2.0 * np.log(n)

    signal = simulate_seeded_binseg_signal(n, true_cps, means, sds, seed)
    map_cps = fit_seeded_binseg(signal)

    print("n =", n)
    print(f"penalty = {penalty:.6f}")
    print("true changepoints =", true_cps)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

"""
Simulate a univariate Gaussian series and detect changepoints with skchange.CROPS.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from skchange.change_detectors import CROPS
from skchange.costs import L2Cost


def simulate_crops_signal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.normal(loc=means[i], scale=sds[i], size=end - start)
    return signal


def fit_crops(signal: np.ndarray) -> tuple[list[int], float, int]:
    detector = CROPS(
        cost=L2Cost(),
        min_penalty=5.0,
        max_penalty=60.0,
        selection_method="bic",
        min_segment_length=10,
        step_size=1,
        split_cost=0.0,
        prune=True,
        pruning_margin=0.0,
    )
    result = detector.fit_predict(pd.DataFrame({"x": signal}))
    return result["ilocs"].to_list(), float(detector.optimal_penalty), len(detector.change_points_lookup)


def main() -> None:
    t0 = time.perf_counter()
    seed = 71
    n = 4000
    true_cps = [0, 800, 1600, 2400, 3200]
    means = [0.0, 4.0, -2.0, 3.5, -1.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]

    signal = simulate_crops_signal(n, true_cps, means, sds, seed)
    map_cps, opt_penalty, n_path = fit_crops(signal)

    print("n =", n)
    print("true changepoints =", true_cps)
    print(f"optimal penalty = {opt_penalty:.6f}")
    print("path solutions =", n_path)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

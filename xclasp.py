"""
Simulate a univariate series and segment it with claspy.BinaryClaSPSegmentation from the local source tree.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

CLASPY_REPO = Path(r"c:/python/public_domain/github/claspy")
if str(CLASPY_REPO) not in sys.path:
    sys.path.insert(0, str(CLASPY_REPO))

from claspy.segmentation import BinaryClaSPSegmentation


def simulate_clasp_signal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise-constant mean.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.normal(loc=means[i], scale=sds[i], size=end - start)
    return signal


def fit_clasp(
    signal: np.ndarray,
    n_segments: int,
    window_size: int,
    k_neighbours: int,
    excl_radius: int,
) -> list[int]:
    """
    Run fixed-K BinaryClaSPSegmentation with a single ClaSP estimator.
    """
    detector = BinaryClaSPSegmentation(
        n_segments=n_segments,
        n_estimators=1,
        window_size=window_size,
        k_neighbours=k_neighbours,
        distance="znormed_euclidean_distance",
        score="roc_auc",
        early_stopping=True,
        validation=None,
        threshold="default",
        excl_radius=excl_radius,
        n_jobs=1,
        random_state=2357,
    )
    return detector.fit_predict(signal).tolist()


def main() -> None:
    t0 = time.perf_counter()
    seed = 101
    n = 500
    n_segments = 5
    window_size = 10
    k_neighbours = 3
    excl_radius = 5
    true_cps = [0, 100, 200, 300, 400]
    means = [0.0, 4.0, -3.0, 3.5, -2.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]

    signal = simulate_clasp_signal(n, true_cps, means, sds, seed)
    estimated = fit_clasp(signal, n_segments, window_size, k_neighbours, excl_radius)

    print("n =", n)
    print("n_segments =", n_segments)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", estimated)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

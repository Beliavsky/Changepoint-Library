"""
Simulate a univariate series and score a single ClaSS split from the local claspy source tree.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

CLASPY_REPO = Path(r"c:/python/public_domain/github/claspy")
if str(CLASPY_REPO) not in sys.path:
    sys.path.insert(0, str(CLASPY_REPO))

from claspy.streaming.clasp import ClaSS


def simulate_class_signal(
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


def fit_class(
    signal: np.ndarray,
    window_size: int,
    k_neighbours: int,
    excl_radius: int,
    threshold: float,
) -> tuple[int | None, float]:
    """
    Run a single-score ClaSS split with score-threshold validation.
    """
    detector = ClaSS(
        window_size=window_size,
        k_neighbours=k_neighbours,
        distance="znormed_euclidean_distance",
        score="f1",
        excl_radius=excl_radius,
    )
    profile = detector.fit_transform(signal)
    cp = detector.split(validation="score_threshold", threshold=threshold)
    score = float(profile[cp]) if cp is not None else float(np.nanmax(profile[np.isfinite(profile)]))
    return cp, score


def main() -> None:
    t0 = time.perf_counter()
    seed = 101
    n = 240
    true_cp = 120
    window_size = 12
    k_neighbours = 3
    excl_radius = 5
    threshold = 0.50
    means = (0.0, 4.0)
    sds = (1.0, 1.0)

    signal = simulate_class_signal(n, true_cp, means, sds, seed)
    estimated, score = fit_class(signal, window_size, k_neighbours, excl_radius, threshold)

    print("n =", n)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("threshold =", threshold)
    print("true changepoint =", true_cp)
    print("estimated changepoint =", estimated)
    print(f"score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

"""
Run BinaryClaSPSegmentation with settings close to changepoynt.CLASP defaults.
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


def simulate_clasp_changepoynt_signal(
    n: int,
    regime_starts: list[int],
    periods: list[int],
    sigma: float,
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate piecewise-sinusoidal signal for CLASP defaults.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        t = np.arange(end - start, dtype=float)
        signal[start:end] = np.sin(2.0 * np.pi * t / periods[i]) + rng.normal(
            loc=0.0, scale=sigma, size=end - start
        )
    return signal


def fit_clasp_changepoynt(signal: np.ndarray) -> BinaryClaSPSegmentation:
    """
    Fit BinaryClaSPSegmentation with the changepoynt wrapper defaults.
    """
    detector = BinaryClaSPSegmentation(
        n_segments="learn",
        n_estimators=10,
        window_size="suss",
        k_neighbours=3,
        distance="znormed_euclidean_distance",
        score="roc_auc",
        early_stopping=True,
        validation="significance_test",
        threshold=1e-15,
        excl_radius=5,
        n_jobs=1,
        random_state=2357,
    )
    detector.fit(signal)
    return detector


def main() -> None:
    t0 = time.perf_counter()
    seed = 123
    n = 880
    true_cps = [0, 220, 440, 660]
    periods = [12, 30, 18, 42]
    sigma = 0.2

    signal = simulate_clasp_changepoynt_signal(n, true_cps, periods, sigma, seed)
    detector = fit_clasp_changepoynt(signal)

    print("n =", n)
    print("n_segments =", "learn")
    print("window_size =", "suss")
    print("n_estimators =", 10)
    print("k_neighbours =", 3)
    print("excl_radius =", 5)
    print("validation =", "significance_test")
    print("threshold =", 1e-15)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", detector.predict().tolist())
    print("window_size used =", detector.window_size)
    print("n_segments used =", detector.n_segments)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

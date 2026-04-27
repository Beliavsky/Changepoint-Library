"""
Simulate a univariate series and score roerich EnergyDistanceCalculator.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROERICH_REPO = Path(r"c:/python/public_domain/github/roerich")
if str(ROERICH_REPO) not in sys.path:
    sys.path.insert(0, str(ROERICH_REPO))

from roerich.change_point import EnergyDistanceCalculator
from roerich.algorithms.cpdc import autoregression_matrix, reference_test


def simulate_roerich_energy_signal(
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


def fit_energy_distance(
    signal: np.ndarray,
    window_size: int,
    periods: int = 1,
    step: int = 1,
) -> tuple[int, float]:
    """
    Run the raw EnergyDistanceCalculator score and return the best shifted changepoint.
    """
    detector = EnergyDistanceCalculator(
        periods=periods,
        window_size=window_size,
        step=step,
        n_runs=1,
    )
    x_auto = autoregression_matrix(signal[:, None], periods=periods, fill_value=0)
    t_score, reference, test = reference_test(x_auto, window_size=window_size, step=1)
    scores = np.array(
        [
            detector.reference_test_predict_n_times(reference[i], test[i])
            for i in range(0, len(reference), step)
        ],
        dtype=float,
    )
    t_eval = np.array([t_score[i] for i in range(0, len(reference), step)], dtype=int)
    best_idx = int(np.argmax(scores))
    return int(t_eval[best_idx] - window_size), float(scores[best_idx])


def main() -> None:
    t0 = time.perf_counter()
    n = 400
    true_cp = 200
    window_size = 40
    periods = 1
    step = 1
    means = (0.0, 3.0)
    sds = (1.0, 1.0)

    signal = simulate_roerich_energy_signal(n, true_cp, means, sds, 101)
    cp, score = fit_energy_distance(signal, window_size, periods, step)

    print("n =", n)
    print("true changepoint =", true_cp)
    print("window_size =", window_size)
    print("periods =", periods)
    print("step =", step)
    print("estimated changepoint =", cp)
    print(f"max raw score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

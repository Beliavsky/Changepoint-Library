"""
Simulate a zero-baseline univariate series and detect collective and point anomalies with skchange.CAPA.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from skchange.anomaly_detectors import CAPA


def simulate_capa_signal(
    n: int,
    segment_anomalies: list[tuple[int, int, float]],
    point_anomalies: list[tuple[int, float]],
    seed: int,
) -> np.ndarray:
    """
    Simulate a zero-baseline Gaussian series with injected segment and point anomalies.
    """
    rng = np.random.default_rng(seed)
    signal = rng.normal(loc=0.0, scale=0.7, size=n)
    for start, end, mean_shift in segment_anomalies:
        signal[start:end] += mean_shift
    for location, value_shift in point_anomalies:
        signal[location] += value_shift
    return signal


def fit_capa(
    signal: np.ndarray,
    segment_penalty: float,
    point_penalty: float,
    min_segment_length: int,
    max_segment_length: int,
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """
    Run CAPA and split the detected anomalies into segment and point intervals.
    """
    detector = CAPA(
        segment_penalty=segment_penalty,
        point_penalty=point_penalty,
        min_segment_length=min_segment_length,
        max_segment_length=max_segment_length,
    )
    result = detector.fit_predict(pd.DataFrame({"x": signal}))
    intervals = [(int(iv.left), int(iv.right)) for iv in result["ilocs"].to_list()]
    segment_intervals = [iv for iv in intervals if iv[1] - iv[0] > 1]
    point_intervals = [iv for iv in intervals if iv[1] - iv[0] == 1]
    return segment_intervals, point_intervals


def main() -> None:
    t0 = time.perf_counter()
    seed = 81
    n = 5_000
    segment_penalty = 24.0
    point_penalty = 14.0
    min_segment_length = 10
    max_segment_length = 500
    true_segments = [(1000, 1100, 5.0), (2800, 2950, -4.5), (4000, 4200, 6.0)]
    true_points = [(2000, 9.0)]

    signal = simulate_capa_signal(n, true_segments, true_points, seed)
    segment_intervals, point_intervals = fit_capa(
        signal,
        segment_penalty,
        point_penalty,
        min_segment_length,
        max_segment_length,
    )

    print("n =", n)
    print(f"segment penalty = {segment_penalty:.6f}")
    print(f"point penalty = {point_penalty:.6f}")
    print("true segment anomalies =", [(start, end) for start, end, _ in true_segments])
    print("true point anomalies =", [(location, location + 1) for location, _ in true_points])
    print("estimated segment anomalies =", segment_intervals)
    print("estimated point anomalies =", point_intervals)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

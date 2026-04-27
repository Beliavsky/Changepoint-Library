"""
Simulate multivariate zero-baseline data and detect subset anomalies with skchange.CAPA.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
from skchange.anomaly_detectors import CAPA


def simulate_mvcapa_signal(
    n: int,
    p: int,
    segment_anomalies: list[tuple[int, int, list[int], float]],
    point_anomalies: list[tuple[int, list[int], float]],
    seed: int,
) -> np.ndarray:
    """
    Simulate zero-baseline Gaussian data with sparse multivariate anomalies.
    """
    rng = np.random.default_rng(seed)
    x = rng.normal(loc=0.0, scale=0.7, size=(n, p))
    for start, end, components, shift in segment_anomalies:
        x[start:end, components] += shift
    for location, components, shift in point_anomalies:
        x[location, components] += shift
    return x


def fit_mvcapa(
    x: np.ndarray,
    segment_penalty: np.ndarray,
    point_penalty: np.ndarray,
    min_segment_length: int,
    max_segment_length: int,
) -> tuple[list[tuple[int, int]], list[list[int]], list[tuple[int, int]], list[list[int]]]:
    """
    Run multivariate CAPA and return segment/point anomalies with affected components.
    """
    detector = CAPA(
        segment_penalty=segment_penalty,
        point_penalty=point_penalty,
        min_segment_length=min_segment_length,
        max_segment_length=max_segment_length,
        find_affected_components=True,
    )
    result = detector.fit_predict(pd.DataFrame(x, columns=[f"x{j + 1}" for j in range(x.shape[1])]))
    intervals = [(int(iv.left), int(iv.right)) for iv in result["ilocs"].to_list()]
    components = [list(map(int, cols)) for cols in result["icolumns"].to_list()]

    segment_intervals = [iv for iv in intervals if iv[1] - iv[0] > 1]
    segment_components = [cols for iv, cols in zip(intervals, components) if iv[1] - iv[0] > 1]
    point_intervals = [iv for iv in intervals if iv[1] - iv[0] == 1]
    point_components = [cols for iv, cols in zip(intervals, components) if iv[1] - iv[0] == 1]
    return segment_intervals, segment_components, point_intervals, point_components


def main() -> None:
    t0 = time.perf_counter()
    seed = 91
    n = 1200
    p = 3
    segment_penalty = np.array([18.0, 26.0, 34.0], dtype=float)
    point_penalty = np.array([14.0, 20.0, 26.0], dtype=float)
    min_segment_length = 5
    max_segment_length = 120
    true_segment_anomalies = [
        (200, 260, [0, 2], 5.0),
        (700, 780, [1], -6.0),
    ]
    true_point_anomalies = [
        (950, [0, 1], 8.0),
    ]

    x = simulate_mvcapa_signal(n, p, true_segment_anomalies, true_point_anomalies, seed)
    seg_iv, seg_cols, pt_iv, pt_cols = fit_mvcapa(
        x,
        segment_penalty,
        point_penalty,
        min_segment_length,
        max_segment_length,
    )

    print("n =", n)
    print("p =", p)
    print("segment penalty =", segment_penalty.tolist())
    print("point penalty =", point_penalty.tolist())
    print("true segment anomalies =", [(s, e, cols) for s, e, cols, _ in true_segment_anomalies])
    print("true point anomalies =", [(t, t + 1, cols) for t, cols, _ in true_point_anomalies])
    print("estimated segment anomalies =", list(zip(seg_iv, seg_cols)))
    print("estimated point anomalies =", list(zip(pt_iv, pt_cols)))
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

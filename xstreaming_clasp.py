"""
Simulate a univariate series and segment it with claspy.StreamingClaSPSegmentation from the local source tree.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

from xclass import simulate_class_signal

CLASPY_REPO = Path(r"c:/python/public_domain/github/claspy")
if str(CLASPY_REPO) not in sys.path:
    sys.path.insert(0, str(CLASPY_REPO))

from claspy.streaming.segmentation import StreamingClaSPSegmentation


def fit_streaming_clasp(
    signal: np.ndarray,
    n_timepoints: int,
    n_warmup: int,
    window_size: int,
    k_neighbours: int,
    jump: int,
    excl_radius: int,
    threshold: float,
) -> tuple[list[int], int]:
    """
    Run StreamingClaSPSegmentation with score-threshold validation and logged changepoints.
    """
    detector = StreamingClaSPSegmentation(
        n_timepoints=n_timepoints,
        n_warmup=n_warmup,
        window_size=window_size,
        k_neighbours=k_neighbours,
        distance="znormed_euclidean_distance",
        score="f1",
        jump=jump,
        validation="score_threshold",
        threshold=threshold,
        log_cps=True,
        excl_radius=excl_radius,
    )
    for value in signal:
        detector.update(float(value))
    return detector.change_points, detector.last_cp


def main() -> None:
    t0 = time.perf_counter()
    seed = 101
    n = 240
    true_cp = 120
    n_timepoints = 200
    n_warmup = 100
    window_size = 12
    k_neighbours = 3
    jump = 5
    excl_radius = 5
    threshold = 0.50

    signal = simulate_class_signal(n, true_cp, (0.0, 4.0), (1.0, 1.0), seed)
    estimated, last_cp = fit_streaming_clasp(
        signal,
        n_timepoints,
        n_warmup,
        window_size,
        k_neighbours,
        jump,
        excl_radius,
        threshold,
    )

    print("n =", n)
    print("n_timepoints =", n_timepoints)
    print("n_warmup =", n_warmup)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("jump =", jump)
    print("excl_radius =", excl_radius)
    print("threshold =", threshold)
    print("true changepoint =", true_cp)
    print("estimated changepoints =", [int(cp) for cp in estimated])
    print("last_cp =", int(last_cp))
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

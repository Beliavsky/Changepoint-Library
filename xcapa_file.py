"""
Read a univariate series from a text file and detect collective and point anomalies with skchange.CAPA.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xcapa import fit_capa


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcapa_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n = signal.shape[0]
    segment_penalty = 24.0
    point_penalty = 14.0
    min_segment_length = 10
    max_segment_length = 500
    segment_intervals, point_intervals = fit_capa(
        signal,
        segment_penalty,
        point_penalty,
        min_segment_length,
        max_segment_length,
    )

    print("file =", data_file)
    print("n =", n)
    print(f"segment penalty = {segment_penalty:.6f}")
    print(f"point penalty = {point_penalty:.6f}")
    print("estimated segment anomalies =", segment_intervals)
    print("estimated point anomalies =", point_intervals)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

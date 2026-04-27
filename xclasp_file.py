"""
Read a univariate series from a file and segment it with claspy.BinaryClaSPSegmentation from local source.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xclasp import fit_clasp


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n_segments = 5
    window_size = 10
    k_neighbours = 3
    excl_radius = 5
    estimated = fit_clasp(signal, n_segments, window_size, k_neighbours, excl_radius)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("n_segments =", n_segments)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("estimated changepoints =", estimated)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

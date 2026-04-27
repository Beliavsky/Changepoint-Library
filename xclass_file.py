"""
Read a univariate series from a file and score a single ClaSS split from local claspy source.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xclass import fit_class


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclass_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    window_size = 12
    k_neighbours = 3
    excl_radius = 5
    threshold = 0.50
    estimated, score = fit_class(signal, window_size, k_neighbours, excl_radius, threshold)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("threshold =", threshold)
    print("estimated changepoint =", estimated)
    print(f"score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

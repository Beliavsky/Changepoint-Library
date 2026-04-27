"""
Read a univariate series from a file and segment it with local-source StreamingClaSPSegmentation.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xstreaming_clasp import fit_streaming_clasp


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xstreaming_clasp_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n_timepoints = 200
    n_warmup = 100
    window_size = 12
    k_neighbours = 3
    jump = 5
    excl_radius = 5
    threshold = 0.50
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

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("n_timepoints =", n_timepoints)
    print("n_warmup =", n_warmup)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("jump =", jump)
    print("excl_radius =", excl_radius)
    print("threshold =", threshold)
    print("estimated changepoints =", [int(cp) for cp in estimated])
    print("last_cp =", int(last_cp))
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

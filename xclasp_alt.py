"""
Second simulated univariate series for claspy.BinaryClaSPSegmentation from the local source tree.
"""

from __future__ import annotations

import time

from xclasp import fit_clasp, simulate_clasp_signal


def main() -> None:
    t0 = time.perf_counter()
    seed = 118
    n = 600
    n_segments = 5
    window_size = 10
    k_neighbours = 3
    excl_radius = 5
    true_cps = [0, 120, 240, 360, 480]
    means = [0.0, 4.0, -3.0, 3.5, -2.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]

    signal = simulate_clasp_signal(n, true_cps, means, sds, seed)
    estimated = fit_clasp(signal, n_segments, window_size, k_neighbours, excl_radius)

    print("n =", n)
    print("n_segments =", n_segments)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", estimated)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

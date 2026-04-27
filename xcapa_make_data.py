"""
Write a zero-baseline univariate anomaly series to a text file for CAPA comparisons.
"""

from __future__ import annotations

import sys

import numpy as np

from xcapa import simulate_capa_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcapa_data.txt"

    seed = 81
    n = 5_000
    true_segments = [(1000, 1100, 5.0), (2800, 2950, -4.5), (4000, 4200, 6.0)]
    true_points = [(2000, 9.0)]
    signal = simulate_capa_signal(n, true_segments, true_points, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# zero-baseline Gaussian anomaly data for skchange.CAPA\n")
        fh.write("# true_segments = " + " ".join(f"[{start},{end})" for start, end, _ in true_segments) + "\n")
        fh.write("# true_points = " + " ".join(f"[{location},{location + 1})" for location, _ in true_points) + "\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true segment anomalies =", [(start, end) for start, end, _ in true_segments])
    print("true point anomalies =", [(location, location + 1) for location, _ in true_points])


if __name__ == "__main__":
    main()

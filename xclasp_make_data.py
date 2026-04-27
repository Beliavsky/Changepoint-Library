"""
Write simulated input data for claspy/BinaryClaSPSegmentation comparisons.
"""

from __future__ import annotations

import sys

from xclasp import simulate_clasp_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_data.txt"

    seed = 101
    n = 500
    true_cps = [0, 100, 200, 300, 400]
    means = [0.0, 4.0, -3.0, 3.5, -2.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]
    signal = simulate_clasp_signal(n, true_cps, means, sds, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# piecewise-Gaussian data for BinaryClaSPSegmentation\n")
        fh.write("# true_cps = " + " ".join(map(str, true_cps)) + "\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

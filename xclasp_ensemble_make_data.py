"""
Write simulated input data for the deterministic ClaSPEnsemble comparison.
"""

from __future__ import annotations

import sys

from xclasp import simulate_clasp_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_ensemble_data.txt"

    seed = 101
    n = 240
    true_cp = 120
    signal = simulate_clasp_signal(n, [0, true_cp], [0.0, 4.0], [1.0, 1.0], seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# piecewise-Gaussian data for deterministic ClaSPEnsemble comparison\n")
        fh.write("# true_cp = 120\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

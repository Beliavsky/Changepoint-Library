"""
Write simulated input data for roerich SlidingWindows comparisons.
"""

from __future__ import annotations

import sys

from xroerich_windows import simulate_roerich_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xroerich_windows_data.txt"

    n = 400
    true_cp = 200
    signal = simulate_roerich_signal(n, true_cp, (0.0, 3.0), (1.0, 1.0), 101)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# piecewise-Gaussian data for roerich SlidingWindows comparison\n")
        fh.write("# true_cp = 200\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

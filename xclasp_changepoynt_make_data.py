"""
Write simulated input data for changepoynt-style CLASP comparisons.
"""

from __future__ import annotations

import sys

from xclasp_changepoynt import simulate_clasp_changepoynt_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_changepoynt_data.txt"

    seed = 123
    n = 880
    true_cps = [0, 220, 440, 660]
    periods = [12, 30, 18, 42]
    sigma = 0.2
    signal = simulate_clasp_changepoynt_signal(n, true_cps, periods, sigma, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# piecewise-sinusoidal data for changepoynt-style CLASP\n")
        fh.write("# true_cps = " + " ".join(map(str, true_cps)) + "\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

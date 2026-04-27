"""
Write simulated input data for changepoints::DP.univar comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_dp_signal(seed: int) -> tuple[np.ndarray, list[int]]:
    """
    Simulate the README-style univariate mean-change signal.
    """
    rng = np.random.default_rng(seed)
    n = 300
    cps = [20, 50, 170]
    means = np.r_[np.repeat(0.0, 20), np.repeat(2.0, 30), np.repeat(0.0, 120), np.repeat(-2.0, 130)]
    y = means + rng.normal(0.0, 1.0, size=n)
    return y, cps


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xdp_univar_data.txt"
    y, cps = simulate_dp_signal(seed=0)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# univariate mean-change data for changepoints::DP.univar\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {len(y)} observations to {data_file}")
    print("true changepoints =", cps)


if __name__ == "__main__":
    main()

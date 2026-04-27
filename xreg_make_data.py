"""
Write simulated regression-change input data for changepoint::cpt.reg(method="AMOC") comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_reg_amoc_data(n: int, true_cp: int, seed: int) -> np.ndarray:
    """
    Simulate a regression series with one changepoint in intercept and slope.
    """
    rng = np.random.default_rng(seed)
    x = np.arange(1, n + 1, dtype=float)
    intercept = np.where(x <= true_cp, 0.0, 50.0)
    slope = np.where(x <= true_cp, 1.0, 0.25)
    y = intercept + slope * x + rng.normal(0.0, 1.0, size=n)
    return np.column_stack([y, np.ones(n, dtype=float), x])


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xreg_data.txt"

    seed = 1
    n = 200
    true_cp = 100
    data = simulate_reg_amoc_data(n, true_cp, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# regression-change data for changepoint::cpt.reg(method='AMOC')\n")
        fh.write(f"# true_cp = {true_cp}\n")
        for row in data:
            fh.write(" ".join(f"{value:.12f}" for value in row) + "\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

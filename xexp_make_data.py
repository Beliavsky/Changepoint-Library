"""
Write simulated input data for changepoint::cpt.meanvar(test.stat="Exponential") comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_exponential_signal(n: int, true_cp: int, mean_left: float, mean_right: float, seed: int) -> np.ndarray:
    """
    Simulate a univariate Exponential series with one changepoint in the mean.
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=float)
    x[:true_cp] = rng.exponential(scale=mean_left, size=true_cp)
    x[true_cp:] = rng.exponential(scale=mean_right, size=n - true_cp)
    return x


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xexp_data.txt"

    seed = 1
    n = 200
    true_cp = 100
    x = simulate_exponential_signal(n, true_cp, 5.0, 12.0, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# single-change Exponential data for changepoint::cpt.meanvar(test.stat='Exponential')\n")
        fh.write(f"# true_cp = {true_cp}\n")
        for value in x:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

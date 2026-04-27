"""
Write simulated input data for changepoint::cpt.meanvar(test.stat="Gamma") comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_gamma_signal(n: int, true_cp: int, shape: float, mean_left: float, mean_right: float, seed: int) -> np.ndarray:
    """
    Simulate a univariate Gamma series with one changepoint in the mean and fixed shape.
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=float)
    x[:true_cp] = rng.gamma(shape=shape, scale=mean_left / shape, size=true_cp)
    x[true_cp:] = rng.gamma(shape=shape, scale=mean_right / shape, size=n - true_cp)
    return x


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xgamma_data.txt"

    seed = 1
    n = 200
    true_cp = 100
    shape = 2.0
    x = simulate_gamma_signal(n, true_cp, shape, 5.0, 12.0, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# single-change Gamma data for changepoint::cpt.meanvar(test.stat='Gamma')\n")
        fh.write(f"# true_cp = {true_cp}\n")
        fh.write(f"# shape = {shape}\n")
        for value in x:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)
    print("shape =", shape)


if __name__ == "__main__":
    main()

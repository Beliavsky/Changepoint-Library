"""
Write simulated input data for changepoint::cpt.meanvar(test.stat="Poisson") comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_poisson_signal(n: int, true_cp: int, rate_left: float, rate_right: float, seed: int) -> np.ndarray:
    """
    Simulate a univariate Poisson series with one changepoint.
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=int)
    x[:true_cp] = rng.poisson(lam=rate_left, size=true_cp)
    x[true_cp:] = rng.poisson(lam=rate_right, size=n - true_cp)
    return x


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xpoisson_data.txt"

    seed = 1
    n = 200
    true_cp = 100
    x = simulate_poisson_signal(n, true_cp, 5.0, 12.0, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# single-change Poisson count data for changepoint::cpt.meanvar(test.stat='Poisson')\n")
        fh.write(f"# true_cp = {true_cp}\n")
        for value in x:
            fh.write(f"{int(value)}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

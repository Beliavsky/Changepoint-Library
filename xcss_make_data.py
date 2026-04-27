"""
Write simulated input data for changepoint::cpt.var(method="AMOC", test.stat="CSS") comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_css_signal(n: int, true_cp: int, sd_left: float, sd_right: float, seed: int) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with one changepoint in variance.
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=float)
    x[:true_cp] = rng.normal(loc=0.0, scale=sd_left, size=true_cp)
    x[true_cp:] = rng.normal(loc=0.0, scale=sd_right, size=n - true_cp)
    return x


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcss_data.txt"

    seed = 1
    n = 200
    true_cp = 100
    x = simulate_css_signal(n, true_cp, 1.0, 10.0, seed)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# single-change Gaussian variance-shift data for changepoint::cpt.var(method='AMOC', test.stat='CSS')\n")
        fh.write(f"# true_cp = {true_cp}\n")
        for value in x:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

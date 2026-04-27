"""
Write simulated input data and thresholds for changepoints::online.univar comparisons.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def simulate_online_signal(seed: int) -> tuple[np.ndarray, list[int]]:
    """
    Simulate a single mean shift for online univariate detection.
    """
    rng = np.random.default_rng(seed)
    n = 150
    cps = [100]
    y = rng.normal(0.0, 1.0, size=n)
    y[100:] += 1.5
    return y, cps


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xonline_univar_data.txt")
    threshold_file = data_file.with_name(data_file.stem + "_thresholds.txt")
    y, cps = simulate_online_signal(seed=0)
    b_vec = np.full(len(y) - 1, 3.5, dtype=float)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# data for changepoints::online.univar with supplied threshold vector\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    with threshold_file.open("w", encoding="utf-8") as fh:
        fh.write("# threshold vector b_t for t >= 2\n")
        for value in b_vec:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {len(y)} observations to {data_file}")
    print(f"wrote {len(b_vec)} thresholds to {threshold_file}")
    print("true changepoints =", cps)


if __name__ == "__main__":
    main()

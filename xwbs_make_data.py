"""
Write simulated input data and WBS intervals for changepoints::WBS.univar comparisons.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def simulate_wbs_signal(seed: int) -> tuple[np.ndarray, list[int]]:
    """
    Simulate the README-style univariate mean-change signal.
    """
    rng = np.random.default_rng(seed)
    n = 300
    cps = [20, 50, 170]
    means = np.r_[np.repeat(0.0, 20), np.repeat(2.0, 30), np.repeat(0.0, 120), np.repeat(-2.0, 130)]
    y = means + rng.normal(0.0, 1.0, size=n)
    return y, cps


def write_intervals(path: Path, n: int, m: int, seed: int) -> None:
    """
    Generate and write WBS intervals as a two-column text file.
    """
    rng = np.random.default_rng(seed)
    alpha = rng.integers(1, n + 1, size=m)
    beta = rng.integers(1, n + 1, size=m)
    lo = np.minimum(alpha, beta)
    hi = np.maximum(alpha, beta)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("# random WBS intervals\n")
        for a, b in zip(lo, hi):
            fh.write(f"{int(a)} {int(b)}\n")


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xwbs_data.txt")
    interval_file = data_file.with_name(data_file.stem + "_intervals.txt")

    seed_signal = 0
    seed_intervals = 1
    m = 300
    y, cps = simulate_wbs_signal(seed_signal)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# univariate mean-change data for changepoints::WBS.univar\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    write_intervals(interval_file, len(y), m, seed_intervals)

    print(f"wrote {len(y)} observations to {data_file}")
    print(f"wrote {m} intervals to {interval_file}")
    print("true changepoints =", cps)


if __name__ == "__main__":
    main()

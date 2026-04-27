"""
Write contaminated input data and WBS intervals for changepoints::WBS.uni.rob comparisons.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from xwbs_make_data import write_intervals


def simulate_wbs_rob_signal(seed: int) -> tuple[np.ndarray, list[int], list[int], list[float]]:
    """
    Simulate a univariate mean-change signal with fixed gross outliers.
    """
    rng = np.random.default_rng(seed)
    n = 300
    cps = [20, 50, 170]
    means = np.r_[np.repeat(0.0, 20), np.repeat(2.0, 30), np.repeat(0.0, 120), np.repeat(-2.0, 130)]
    y = means + rng.normal(0.0, 1.0, size=n)
    outlier_idx = [15, 35, 60, 90, 140, 180, 220, 260, 290]
    outlier_shift = [15.0, -12.0, 14.0, -16.0, 13.0, 12.0, -14.0, 15.0, -13.0]
    y[np.array(outlier_idx) - 1] += np.array(outlier_shift)
    return y, cps, outlier_idx, outlier_shift


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xwbs_rob_data.txt")
    interval_file = data_file.with_name(data_file.stem + "_intervals.txt")

    seed_signal = 0
    seed_intervals = 1
    m = 300
    y, cps, outlier_idx, outlier_shift = simulate_wbs_rob_signal(seed_signal)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# contaminated univariate mean-change data for changepoints::WBS.uni.rob\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        fh.write("# outlier_idx = " + " ".join(str(i) for i in outlier_idx) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    write_intervals(interval_file, len(y), m, seed_intervals)

    print(f"wrote {len(y)} observations to {data_file}")
    print(f"wrote {m} intervals to {interval_file}")
    print("true changepoints =", cps)
    print("outlier_idx =", outlier_idx)
    print("outlier_shift =", outlier_shift)


if __name__ == "__main__":
    main()

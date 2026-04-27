"""
Write contaminated input data for changepoints::ARC comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_arc_signal(seed: int) -> tuple[np.ndarray, list[int], list[int], list[float]]:
    """
    Simulate a single-change Gaussian mean-shift signal with fixed gross outliers.
    """
    rng = np.random.default_rng(seed)
    n = 1000
    cps = [500]
    y = np.r_[rng.normal(0.0, 1.0, 500), rng.normal(1.0, 1.0, 500)]
    outlier_idx = [80, 150, 220, 310, 420, 560, 640, 730, 820, 910]
    outlier_shift = [12.0, -11.0, 10.0, -13.0, 11.0, 12.0, -10.0, 13.0, -12.0, 11.0]
    y[np.array(outlier_idx) - 1] += np.array(outlier_shift)
    return y, cps, outlier_idx, outlier_shift


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xarc_data.txt"
    y, cps, outlier_idx, outlier_shift = simulate_arc_signal(seed=0)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# contaminated univariate mean-change data for changepoints::ARC\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        fh.write("# outlier_idx = " + " ".join(str(i) for i in outlier_idx) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {len(y)} observations to {data_file}")
    print("true changepoints =", cps)
    print("outlier_idx =", outlier_idx)
    print("outlier_shift =", outlier_shift)


if __name__ == "__main__":
    main()

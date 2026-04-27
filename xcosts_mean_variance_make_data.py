"""
Create a univariate data file with piecewise changes in both mean and variance.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_mean_variance_signal(
    n: int,
    bkps: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise mean and variance.
    """
    if len(bkps) != len(means) or len(bkps) != len(sds):
        raise ValueError("bkps, means, and sds must have the same length")
    if bkps[-1] != n:
        raise ValueError("final breakpoint must equal n")

    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    start = 0
    for end, mean, sd in zip(bkps, means, sds):
        signal[start:end] = rng.normal(loc=mean, scale=sd, size=end - start)
        start = end
    return signal


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xcosts_mean_variance_data.txt"
    seed = 11
    n = 600
    true_bkps = [150, 300, 450, n]
    means = [0.0, 2.5, -1.0, 1.5]
    sds = [1.0, 2.2, 0.6, 1.8]

    signal = simulate_mean_variance_signal(n, true_bkps, means, sds, seed)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        f.write("# means = " + " ".join(f"{x:.6f}" for x in means) + "\n")
        f.write("# sds = " + " ".join(f"{x:.6f}" for x in sds) + "\n")
        np.savetxt(f, signal, fmt="%.6f")

    print(f"wrote {n} observations to {outfile}")
    print("true bkps =", true_bkps)


if __name__ == "__main__":
    main()

"""
Create a univariate Gaussian data file for changepoint.ArgpCpd comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_argpcpd_signal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise mean and variance.
    """
    if len(regime_starts) != len(means) or len(regime_starts) != len(sds):
        raise ValueError("regime_starts, means, and sds must have the same length")
    if regime_starts[0] != 0:
        raise ValueError("regime_starts must begin at 0")

    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.normal(loc=means[i], scale=sds[i], size=end - start)
    return signal


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xargpcpd_data.txt"
    seed = 31
    n = 600
    true_bkps = [0, 200, 400]
    means = [0.0, 4.0, -3.0]
    sds = [0.4, 0.4, 0.5]

    signal = simulate_argpcpd_signal(n, true_bkps, means, sds, seed)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        f.write("# means = " + " ".join(f"{x:.6f}" for x in means) + "\n")
        f.write("# sds = " + " ".join(f"{x:.6f}" for x in sds) + "\n")
        np.savetxt(f, signal, fmt="%.6f")

    print(f"wrote {n} observations to {outfile}")
    print("true changepoints =", true_bkps)


if __name__ == "__main__":
    main()

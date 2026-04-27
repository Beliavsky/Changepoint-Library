"""
Create a piecewise-Poisson count data file for changepoint.Bocpd comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_poisson_signal(
    n: int,
    regime_starts: list[int],
    rates: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a piecewise-Poisson count series.
    """
    if len(regime_starts) != len(rates):
        raise ValueError("regime_starts and rates must have the same length")
    if regime_starts[0] != 0:
        raise ValueError("regime_starts must begin at 0")

    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=np.int64)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.poisson(lam=rates[i], size=end - start)
    return signal


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xbocpd_poisson_data.txt"
    seed = 41
    n = 600
    true_bkps = [0, 200, 400]
    rates = [2.0, 9.0, 4.0]

    signal = simulate_poisson_signal(n, true_bkps, rates, seed)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        f.write("# rates = " + " ".join(f"{x:.6f}" for x in rates) + "\n")
        np.savetxt(f, signal, fmt="%d")

    print(f"wrote {n} observations to {outfile}")
    print("true changepoints =", true_bkps)


if __name__ == "__main__":
    main()

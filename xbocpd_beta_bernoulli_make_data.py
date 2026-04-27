"""
Create a piecewise-Bernoulli binary data file for changepoint.Bocpd comparisons.
"""

from __future__ import annotations

import sys

import numpy as np


def simulate_bernoulli_signal(
    n: int,
    regime_starts: list[int],
    probs: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a piecewise-Bernoulli binary series.
    """
    if len(regime_starts) != len(probs):
        raise ValueError("regime_starts and probs must have the same length")
    if regime_starts[0] != 0:
        raise ValueError("regime_starts must begin at 0")

    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=np.int64)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.binomial(n=1, p=probs[i], size=end - start)
    return signal


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xbocpd_beta_bernoulli_data.txt"
    seed = 51
    n = 600
    true_bkps = [0, 200, 400]
    probs = [0.15, 0.80, 0.35]

    signal = simulate_bernoulli_signal(n, true_bkps, probs, seed)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        f.write("# probs = " + " ".join(f"{x:.6f}" for x in probs) + "\n")
        np.savetxt(f, signal, fmt="%d")

    print(f"wrote {n} observations to {outfile}")
    print("true changepoints =", true_bkps)


if __name__ == "__main__":
    main()

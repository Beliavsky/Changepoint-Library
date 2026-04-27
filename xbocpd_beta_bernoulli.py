"""
Simulate a piecewise-Bernoulli binary series and detect changepoints with changepoint.Bocpd.
"""

from __future__ import annotations

import time

import changepoint as cpt
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
    t0 = time.perf_counter()
    seed = 51
    n = 600
    true_cps = [0, 200, 400]
    probs = [0.15, 0.80, 0.35]

    signal = simulate_bernoulli_signal(n, true_cps, probs, seed)
    bocpd = cpt.Bocpd(cpt.BetaBernoulli(), 120.0)
    rs = [bocpd.step(bool(x)) for x in signal]
    map_cps = cpt.map_changepoints(rs)

    print("n =", n)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

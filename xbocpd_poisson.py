"""
Simulate a piecewise-Poisson count series and detect changepoints with changepoint.Bocpd.
"""

from __future__ import annotations

import time

import changepoint as cpt
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
    t0 = time.perf_counter()
    seed = 41
    n = 600
    true_cps = [0, 200, 400]
    rates = [2.0, 9.0, 4.0]

    signal = simulate_poisson_signal(n, true_cps, rates, seed)
    bocpd = cpt.Bocpd(cpt.PoissonGamma(), 120.0)
    rs = [bocpd.step(int(x)) for x in signal]
    map_cps = cpt.map_changepoints(rs)

    print("n =", n)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

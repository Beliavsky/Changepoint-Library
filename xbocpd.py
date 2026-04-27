"""
Simulate a univariate Gaussian series and detect changepoints with changepoint.Bocpd.
"""

from __future__ import annotations

import time

import changepoint as cpt
import numpy as np


def simulate_bocpd_signal(
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
    t0 = time.perf_counter()
    seed = 21
    n = 600
    lam = 120.0
    true_cps = [0, 200, 400]
    means = [0.0, 3.5, -1.5]
    sds = [1.0, 1.0, 1.6]

    signal = simulate_bocpd_signal(n, true_cps, means, sds, seed)
    bocpd = cpt.Bocpd(cpt.NormalGamma(), lam)
    rs = [bocpd.step(float(x)) for x in signal]
    map_cps = cpt.map_changepoints(rs)

    print("n =", n)
    print("lam =", lam)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

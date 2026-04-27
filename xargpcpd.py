"""
Simulate a univariate Gaussian series and detect changepoints with changepoint.ArgpCpd.
"""

from __future__ import annotations

import time

import changepoint as cpt
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


def fit_argpcpd(signal: np.ndarray) -> list[int]:
    """
    Run ArgpCpd and return MAP changepoints.
    """
    argp = cpt.ArgpCpd(
        logistic_hazard_h=-5.0,
        scale=3.0,
        noise_level=0.01,
        max_lag=12,
    )
    rs = [argp.step(float(x)) for x in signal]
    return cpt.map_changepoints(rs)


def main() -> None:
    t0 = time.perf_counter()
    seed = 31
    n = 600
    true_cps = [0, 200, 400]
    means = [0.0, 4.0, -3.0]
    sds = [0.4, 0.4, 0.5]

    signal = simulate_argpcpd_signal(n, true_cps, means, sds, seed)
    map_cps = fit_argpcpd(signal)

    print("n =", n)
    print("true changepoints =", true_cps)
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

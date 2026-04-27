#!/usr/bin/env python
"""Standalone demo for the local `bocd` package."""

from __future__ import annotations

import time

from xbocd_utils import run_bocd_student_t, simulate_piecewise_normal


def main() -> None:
    n = 400
    lam = 100.0
    regime_starts = [0, 150, 300]
    means = [0.0, 3.0, -2.0]
    sds = [1.0, 1.0, 1.0]
    seed = 20260422

    x, true_cps = simulate_piecewise_normal(n, regime_starts, means, sds, seed)
    t0 = time.perf_counter()
    out = run_bocd_student_t(x, lam=lam, mu=0.0, kappa=1.0, alpha=1.0, beta=1.0)
    elapsed = time.perf_counter() - t0

    print(f"n = {n}")
    print(f"lambda = {lam}")
    print(f"true changepoints = {true_cps}")
    print(f"estimated changepoints = {out['map_cps']}")
    print(f"final MAP run length = {out['final_rt']}")
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

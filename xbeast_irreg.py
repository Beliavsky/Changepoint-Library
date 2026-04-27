#!/usr/bin/env python
"""Standalone irregular trend-only BEAST-style demo."""

from __future__ import annotations

import time

from xbeast_utils import (
    simulate_irregular_time_grid,
    simulate_piecewise_linear_irregular,
    solve_beast_irreg_trend_only,
    top_cp_probabilities,
)


def main() -> None:
    n = 220
    max_cp = 4
    min_seg_len = 28
    regime_starts = [0, 70, 145]
    intercepts = [1.0, 6.0, -2.0]
    slopes = [0.05, -0.01, 0.03]
    sigma = 0.30
    seed_time = 20260423
    seed_y = 20260424

    t = simulate_irregular_time_grid(n, start_time=0.0, dt_mean=1.0, dt_jitter=0.35, seed=seed_time)
    y, true_cps = simulate_piecewise_linear_irregular(t, regime_starts, intercepts, slopes, sigma, seed_y)

    t0 = time.perf_counter()
    out = solve_beast_irreg_trend_only(t, y, max_cp=max_cp, min_seg_len=min_seg_len)
    elapsed = time.perf_counter() - t0

    print(f"n = {n}")
    print("season = none")
    print("irregular = True")
    print(f"max_cp = {max_cp}")
    print(f"min_seg_len = {min_seg_len}")
    print(f"true changepoints = {true_cps}")
    print(f"estimated changepoints = {out['best_cps']}")
    print("model weights =", [round(float(w), 6) for w in out["model_weights"]])
    print("top cp probabilities =", [(i, round(p, 6)) for i, p in top_cp_probabilities(out["cp_prob"])])
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Standalone regular trend-only BEAST-style demo with season='none'."""

from __future__ import annotations

import time

from xbeast_utils import simulate_piecewise_linear, solve_beast_trend_only, top_cp_probabilities


def main() -> None:
    n = 240
    max_cp = 4
    min_seg_len = 30
    regime_starts = [0, 80, 160]
    intercepts = [1.0, 4.6, -1.2]
    slopes = [0.03, -0.02, 0.015]
    sigma = 0.35
    seed = 20260421

    y, true_cps = simulate_piecewise_linear(n, regime_starts, intercepts, slopes, sigma, seed)
    t0 = time.perf_counter()
    out = solve_beast_trend_only(y, max_cp=max_cp, min_seg_len=min_seg_len)
    elapsed = time.perf_counter() - t0

    print(f"n = {n}")
    print("season = none")
    print(f"max_cp = {max_cp}")
    print(f"min_seg_len = {min_seg_len}")
    print(f"true changepoints = {true_cps}")
    print(f"estimated changepoints = {out['best_cps']}")
    print("model weights =", [round(float(w), 6) for w in out["model_weights"]])
    print("top cp probabilities =", [(i, round(p, 6)) for i, p in top_cp_probabilities(out["cp_prob"])])
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

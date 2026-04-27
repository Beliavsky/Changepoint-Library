#!/usr/bin/env python
"""Standalone regular multi-series BEAST123-style demo."""

from __future__ import annotations

import time

import numpy as np

from xbeast_utils import simulate_piecewise_linear, solve_beast123_regular, top_cp_probabilities


def main() -> None:
    n = 240
    p = 3
    max_cp = 4
    min_seg_len = 30
    regime_starts = [0, 80, 160]
    sigma = 0.35
    intercepts = np.array(
        [
            [1.0, -1.0, 2.5],
            [4.6, 2.0, -0.5],
            [-1.2, 3.2, 1.0],
        ],
        dtype=float,
    )
    slopes = np.array(
        [
            [0.03, 0.015, -0.02],
            [-0.02, -0.010, 0.025],
            [0.015, 0.020, -0.015],
        ],
        dtype=float,
    )
    x = np.empty((n, p), dtype=float)
    true_cps: list[int] = []
    for j in range(p):
        x[:, j], true_cps = simulate_piecewise_linear(
            n, regime_starts, intercepts[:, j].tolist(), slopes[:, j].tolist(), sigma, 20260430 + j
        )

    t0 = time.perf_counter()
    out = solve_beast123_regular(x, max_cp=max_cp, min_seg_len=min_seg_len)
    elapsed = time.perf_counter() - t0

    print(f"n = {n}")
    print(f"p = {p}")
    print("season = none")
    print("beast123 regular interface")
    print(f"max_cp = {max_cp}")
    print(f"min_seg_len = {min_seg_len}")
    print(f"true changepoints = {true_cps}")
    for j, cps in enumerate(out["best_cps"], start=1):
        print(f"series {j} changepoints = {cps}")
    print("top mean cp probabilities =", [(i, round(pv, 6)) for i, pv in top_cp_probabilities(out["mean_cp_prob"])])
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

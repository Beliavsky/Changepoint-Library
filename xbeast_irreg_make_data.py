#!/usr/bin/env python
"""Write an irregular trend-only BEAST-style dataset to a text file."""

from __future__ import annotations

import sys

from xbeast_utils import (
    simulate_irregular_time_grid,
    simulate_piecewise_linear_irregular,
    write_time_series_file,
)


def main() -> None:
    output_path = sys.argv[1] if len(sys.argv) > 1 else "xbeast_irreg_data.txt"
    n = 220
    regime_starts = [0, 70, 145]
    intercepts = [1.0, 6.0, -2.0]
    slopes = [0.05, -0.01, 0.03]
    sigma = 0.30
    seed_time = 20260423
    seed_y = 20260424

    t = simulate_irregular_time_grid(n, start_time=0.0, dt_mean=1.0, dt_jitter=0.35, seed=seed_time)
    y, true_cps = simulate_piecewise_linear_irregular(t, regime_starts, intercepts, slopes, sigma, seed_y)
    write_time_series_file(output_path, t, y, true_cps)
    print(f"wrote {n} observations to {output_path}")
    print(f"true changepoints = {true_cps}")


if __name__ == "__main__":
    main()

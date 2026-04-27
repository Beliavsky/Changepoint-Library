#!/usr/bin/env python
"""Write a regular harmonic-seasonal BEAST-style dataset to a text file."""

from __future__ import annotations

import sys

from xbeast_utils import simulate_piecewise_linear_harmonic, write_series_file


def main() -> None:
    output_path = sys.argv[1] if len(sys.argv) > 1 else "xbeast_seasonal_data.txt"
    n = 240
    period = 12.0
    regime_starts = [0, 80, 160]
    intercepts = [1.0, 4.6, -1.2]
    slopes = [0.03, -0.02, 0.015]
    sigma = 0.35
    amp_sin = 1.4
    amp_cos = -0.8
    seed = 20260422

    y, true_cps = simulate_piecewise_linear_harmonic(
        n, regime_starts, intercepts, slopes, sigma, period, amp_sin, amp_cos, seed
    )
    write_series_file(output_path, y, true_cps)
    print(f"wrote {n} observations to {output_path}")
    print(f"true changepoints = {true_cps}")


if __name__ == "__main__":
    main()

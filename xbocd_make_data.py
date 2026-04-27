#!/usr/bin/env python
"""Write a one-column data file for the local `bocd` package workflow."""

from __future__ import annotations

import sys

from xbocd_utils import simulate_piecewise_normal, write_series_file


def main() -> None:
    output_path = sys.argv[1] if len(sys.argv) > 1 else "xbocd_data.txt"
    n = 400
    regime_starts = [0, 150, 300]
    means = [0.0, 3.0, -2.0]
    sds = [1.0, 1.0, 1.0]
    seed = 20260422

    x, true_cps = simulate_piecewise_normal(n, regime_starts, means, sds, seed)
    write_series_file(output_path, x, true_cps)
    print(f"wrote {n} observations to {output_path}")
    print(f"true changepoints = {true_cps}")


if __name__ == "__main__":
    main()

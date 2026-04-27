#!/usr/bin/env python
"""Write a regular multi-series BEAST123-style dataset to a text file."""

from __future__ import annotations

import sys

import numpy as np

from xbeast_utils import simulate_piecewise_linear, write_matrix_file


def main() -> None:
    output_path = sys.argv[1] if len(sys.argv) > 1 else "xbeast123_data.txt"
    n = 240
    p = 3
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

    write_matrix_file(output_path, x, true_cps)
    print(f"wrote {n}x{p} observations to {output_path}")
    print(f"true changepoints = {true_cps}")


if __name__ == "__main__":
    main()

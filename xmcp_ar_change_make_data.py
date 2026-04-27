#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_ar_change_data.txt")
    rng = np.random.default_rng(42)

    n = 180
    x = np.arange(1, n + 1, dtype=float)
    mu = np.empty(n, dtype=float)
    mu[:80] = -20.0
    mu[80:140] = -20.0 + 0.5 * (x[80:140] - 80.0)
    mu[140:] = -20.0 + 0.5 * (140.0 - 80.0)
    sigma = 5.0

    y = np.empty(n, dtype=float)
    eps = rng.normal(0.0, sigma, size=n)
    y[0] = mu[0] + eps[0]
    y[1] = mu[1] + 0.7 * (y[0] - mu[0]) + eps[1]
    for i in range(2, 80):
        y[i] = mu[i] + 0.7 * (y[i - 1] - mu[i - 1]) - 0.4 * (y[i - 2] - mu[i - 2]) + eps[i]
    for i in range(80, n):
        ar1 = 0.5 - 0.02 * (x[i] - 80.0)
        y[i] = mu[i] + ar1 * (y[i - 1] - mu[i - 1]) + eps[i]

    data = np.column_stack([x, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {n} rows")
    print("true params = {cp_1: 80, cp_2: 140, int_1: -20, sigma_1: 5, ar1_1: 0.7, ar2_1: -0.4, x_2: 0.5, ar1_2: 0.5, ar1_x_2: -0.02}")


if __name__ == "__main__":
    main()

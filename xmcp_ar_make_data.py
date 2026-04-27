#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_ar_data.txt")
    rng = np.random.default_rng(42)

    n = 120
    time = np.arange(1, n + 1, dtype=float)
    phi = 0.7
    sigma = 5.0
    mu = 20.0
    eps = rng.normal(0.0, sigma, size=n)
    y = np.empty(n, dtype=float)
    y[0] = mu + eps[0]
    for i in range(1, n):
        y[i] = mu + phi * (y[i - 1] - mu) + eps[i]

    data = np.column_stack([time, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {n} rows")
    print("true params = {int_1: 20, ar1_1: 0.7, sigma_1: 5}")


if __name__ == "__main__":
    main()

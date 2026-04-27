#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_arsigma_data.txt")
    rng = np.random.default_rng(123)

    n = 140
    x = np.arange(1, n + 1, dtype=float)
    mu = np.empty(n, dtype=float)
    mu[:70] = 15.0
    mu[70:] = 15.0 + 0.3 * (x[70:] - 70.0)

    sigma = np.empty(n, dtype=float)
    sigma[:70] = np.maximum(0.0, 3.0 + 0.03 * x[:70])
    sigma[70:] = 6.0

    z = rng.normal(size=n)
    eps = np.zeros(n, dtype=float)
    eps[0] = sigma[0] * z[0]
    for i in range(1, 70):
        eps[i] = 0.5 * eps[i - 1] + sigma[i] * z[i]
    eps[70] = sigma[70] * z[70]
    if n > 71:
        eps[71] = (-0.2 + 0.003 * x[71]) * eps[70] + 0.15 * eps[69] + sigma[71] * z[71]
    for i in range(72, n):
        eps[i] = (-0.2 + 0.003 * x[i]) * eps[i - 1] + 0.15 * eps[i - 2] + sigma[i] * z[i]

    y = mu + eps
    data = np.column_stack([x, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {n} rows")
    print("true params = {cp_1: 70, int_1: 15, x_2: 0.3, sigma_1: 3, sigma_x_1: 0.03, sigma_2: 6, ar1_1: 0.5, ar1_2: -0.2, ar1_x_2: 0.003, ar2_2: 0.15, ar2_x_2: 0.0}")


if __name__ == "__main__":
    main()

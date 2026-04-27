#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_fit_data.txt")
    rng = np.random.default_rng(123)

    n = 100
    time = np.arange(1, n + 1, dtype=float)
    mean = np.empty(n, dtype=float)
    mean[:30] = 10.0
    mean[30:65] = 10.0 + 0.45 * (time[30:65] - 30.0)
    mean[65:] = 18.0 + 0.35 * (time[65:] - 65.0)
    response = mean + rng.normal(0.0, 4.0, size=n)

    data = np.column_stack([time, response])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {n} rows")
    print("true changepoints = [30, 65]")


if __name__ == "__main__":
    main()

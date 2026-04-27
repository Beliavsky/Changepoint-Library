#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xcpm_exponential_adjusted_data.txt")
    n = 180
    cps = [60, 120]
    rng = np.random.default_rng(2418)
    x = np.empty(n)
    x[:60] = rng.exponential(scale=0.6, size=60)
    x[60:120] = rng.exponential(scale=1.6, size=60)
    x[120:] = rng.exponential(scale=0.8, size=60)
    np.savetxt(out_path, x, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with {n} rows")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

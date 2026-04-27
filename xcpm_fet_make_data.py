#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xcpm_fet_data.txt")
    n = 180
    cps = [60, 120]
    rng = np.random.default_rng(2409)
    x = np.empty(n, dtype=int)
    x[:60] = rng.binomial(1, 0.15, size=60)
    x[60:120] = rng.binomial(1, 0.65, size=60)
    x[120:] = rng.binomial(1, 0.30, size=60)
    np.savetxt(out_path, x, fmt="%d", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with {n} rows")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

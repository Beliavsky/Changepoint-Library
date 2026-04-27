#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xchangepointnp_data.txt")
    n = 120
    cps = [40, 80]
    rng = np.random.default_rng(123)
    x = np.empty(n)
    x[:40] = rng.normal(0.0, 1.0, size=40)
    x[40:80] = rng.normal(2.5, 1.0, size=40)
    x[80:] = rng.normal(-1.0, 1.0, size=40)
    np.savetxt(out_path, x, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with {n} rows")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xstrucchange_fstats_data.txt")
    n = 120
    cps = [40, 80]
    rng = np.random.default_rng(913)
    x = np.linspace(-1.0, 1.0, n)
    beta0 = np.empty(n)
    beta1 = np.empty(n)
    beta0[:40] = -0.5
    beta1[:40] = 1.2
    beta0[40:80] = 1.0
    beta1[40:80] = -0.8
    beta0[80:] = -1.2
    beta1[80:] = 0.9
    y = beta0 + beta1 * x + rng.normal(0.0, 0.18, size=n)
    dat = np.column_stack([y, x])
    np.savetxt(out_path, dat, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with {n} rows and {dat.shape[1]} columns")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

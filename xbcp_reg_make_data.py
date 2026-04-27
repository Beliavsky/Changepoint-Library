#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xbcp_reg_data.txt")
    n = 120
    cps = [40, 80]
    x = np.linspace(-1.0, 1.0, n)
    a = np.concatenate([np.full(40, 0.5), np.full(40, 2.0), np.full(40, -1.0)])
    b = np.concatenate([np.full(40, 1.5), np.full(40, -2.0), np.full(40, 0.8)])
    rng = np.random.default_rng(123)
    y = a + b * x + rng.normal(0.0, 0.25, size=n)
    out = np.column_stack([y, x])
    np.savetxt(out_path, out, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with shape {out.shape[0]} x {out.shape[1]}")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xstrucchange_mefp_me_data.txt")
    n = 120
    history_n = 50
    cp = 78
    rng = np.random.default_rng(1409)
    x = np.linspace(-1.0, 1.0, n)
    beta0 = np.full(n, -0.5)
    beta1 = np.full(n, 1.1)
    beta0[cp:] = 0.9
    beta1[cp:] = -0.7
    y = beta0 + beta1 * x + rng.normal(0.0, 0.18, size=n)
    dat = np.column_stack([y, x])
    np.savetxt(out_path, dat, fmt="%.12f", header=f"history_n = {history_n} true_bkps = {cp}")
    print(f"wrote {out_path.name} with {n} rows and {dat.shape[1]} columns")
    print(f"history_n = {history_n}")
    print(f"true changepoints = [{cp}]")


if __name__ == "__main__":
    main()

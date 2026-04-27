#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xecp_cp3o_delta_data.txt")
    n1 = n2 = n3 = 40
    cps = [40, 80]
    rng = np.random.default_rng(2601)
    x1 = rng.normal(loc=(0.0, 0.0), scale=(0.25, 0.35), size=(n1, 2))
    x2 = rng.normal(loc=(3.0, 3.0), scale=(0.30, 0.25), size=(n2, 2))
    x3 = rng.normal(loc=(-2.0, 1.5), scale=(0.20, 0.30), size=(n3, 2))
    x = np.vstack([x1, x2, x3])
    np.savetxt(out_path, x, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with shape {x.shape[0]} x {x.shape[1]}")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

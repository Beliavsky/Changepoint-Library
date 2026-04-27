#!/usr/bin/env python
"""Generate a multivariate series for changepoints::BS.cov."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) >= 2 else Path("xbs_cov_data.txt")
    seed = 1234
    n = 180
    p = 4
    true_bkps = [60, 120]

    cov1 = np.array(
        [
            [1.0, 0.7, 0.0, 0.0],
            [0.7, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.2],
            [0.0, 0.0, 0.2, 1.0],
        ]
    )
    cov2 = np.array(
        [
            [1.8, 0.1, 0.0, 0.0],
            [0.1, 0.7, 0.0, 0.0],
            [0.0, 0.0, 0.8, -0.5],
            [0.0, 0.0, -0.5, 1.5],
        ]
    )
    cov3 = np.array(
        [
            [0.9, -0.6, 0.2, 0.0],
            [-0.6, 1.4, 0.0, 0.0],
            [0.2, 0.0, 1.2, 0.6],
            [0.0, 0.0, 0.6, 1.1],
        ]
    )

    rng = np.random.default_rng(seed)
    x1 = rng.multivariate_normal(np.zeros(p), cov1, size=true_bkps[0])
    x2 = rng.multivariate_normal(np.zeros(p), cov2, size=true_bkps[1] - true_bkps[0])
    x3 = rng.multivariate_normal(np.zeros(p), cov3, size=n - true_bkps[1])
    x = np.vstack((x1, x2, x3))

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write(f"# seed = {seed}\n")
        fh.write(f"# n = {n}\n")
        fh.write(f"# p = {p}\n")
        fh.write("# true_bkps = " + " ".join(str(v) for v in true_bkps) + "\n")
        fh.write("# rows are time points and columns are variables\n")
        np.savetxt(fh, x, fmt="%.12f")

    print(f"wrote {data_file} with shape {x.shape[0]} x {x.shape[1]}")


if __name__ == "__main__":
    main()

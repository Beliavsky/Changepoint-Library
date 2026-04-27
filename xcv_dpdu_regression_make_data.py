#!/usr/bin/env python
"""Write data for changepoints::CV.search.DPDU.regression comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xcv_dpdu_regression_data.txt")
    rng = np.random.default_rng(123)

    n = 120
    p = 3
    true_cps = [40, 80]
    x = rng.normal(size=(n, p))
    beta1 = np.array([1.5, 0.0, -1.0])
    beta2 = np.array([0.2, 1.8, -0.2])
    beta3 = np.array([-1.2, 0.5, 1.1])
    y = np.empty(n)
    y[:40] = x[:40] @ beta1 + rng.normal(scale=0.5, size=40)
    y[40:80] = x[40:80] @ beta2 + rng.normal(scale=0.5, size=40)
    y[80:] = x[80:] @ beta3 + rng.normal(scale=0.5, size=40)

    mat = np.column_stack((y, x))
    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# y x1 x2 x3 for changepoints::CV.search.DPDU.regression\n")
        fh.write(f"# n = {n}\n")
        fh.write(f"# p = {p}\n")
        fh.write("# true_bkps = " + " ".join(str(v) for v in true_cps) + "\n")
        np.savetxt(fh, mat, fmt="%.12f")

    print(f"wrote {data_file} with {n} rows and {p + 1} columns")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

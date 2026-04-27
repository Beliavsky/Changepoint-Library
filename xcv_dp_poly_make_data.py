#!/usr/bin/env python
"""Write data for changepoints::CV.search.DP.poly comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def eval_poly(coef: np.ndarray, t: np.ndarray, n: int) -> np.ndarray:
    x = t / n
    out = np.zeros_like(x, dtype=float)
    power = np.ones_like(x, dtype=float)
    for c in coef:
        out += c * power
        power *= x
    return out


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xcv_dp_poly_data.txt")
    rng = np.random.default_rng(123)

    n = 150
    r = 2
    true_cps = [50, 100]
    t = np.arange(1, n + 1, dtype=float)
    coef1 = np.array([0.0, 1.5, -2.5])
    coef2 = np.array([2.5, -8.0, 10.0])
    coef3 = np.array([-2.0, 5.0, -2.5])
    y = np.empty(n, dtype=float)
    y[:50] = eval_poly(coef1, t[:50], n)
    y[50:100] = eval_poly(coef2, t[50:100], n)
    y[100:] = eval_poly(coef3, t[100:], n)
    y += rng.normal(scale=0.05, size=n)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# y for changepoints::CV.search.DP.poly\n")
        fh.write(f"# n = {n}\n")
        fh.write(f"# r = {r}\n")
        fh.write("# true_bkps = " + " ".join(str(v) for v in true_cps) + "\n")
        np.savetxt(fh, y, fmt="%.12f")

    print(f"wrote {data_file} with {n} rows")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

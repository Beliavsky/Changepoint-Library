#!/usr/bin/env python
"""Write data for changepoints::DP.VAR1 comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xdp_var1_data.txt")
    rng = np.random.default_rng(123)

    p = 2
    n_trans = 120
    true_cps = [40, 80]
    a1 = np.array([[0.6, 0.0], [0.0, 0.4]])
    a2 = np.array([[-0.4, 0.3], [0.2, 0.5]])
    a3 = np.array([[0.7, -0.2], [0.1, -0.5]])
    data = np.zeros((p, n_trans + 1), dtype=float)
    data[:, 0] = np.array([0.2, -0.1])

    for t in range(n_trans):
        if t < 40:
            a = a1
        elif t < 80:
            a = a2
        else:
            a = a3
        data[:, t + 1] = a @ data[:, t] + rng.normal(scale=0.1, size=p)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write("# DATA matrix for changepoints::DP.VAR1 (rows = time, cols = variables)\n")
        fh.write(f"# p = {p}\n")
        fh.write(f"# n_transitions = {n_trans}\n")
        fh.write("# true_bkps = " + " ".join(str(v) for v in true_cps) + "\n")
        np.savetxt(fh, data.T, fmt="%.12f")

    print(f"wrote {data_file} with shape {data.shape[1]} x {data.shape[0]}")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

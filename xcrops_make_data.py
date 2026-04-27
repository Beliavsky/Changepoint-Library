"""
Create a univariate Gaussian data file for skchange.CROPS comparisons.
"""

from __future__ import annotations

import sys

import numpy as np

from xcrops import simulate_crops_signal


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xcrops_data.txt"
    seed = 71
    n = 4000
    true_bkps = [0, 800, 1600, 2400, 3200]
    means = [0.0, 4.0, -2.0, 3.5, -1.5]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0]

    signal = simulate_crops_signal(n, true_bkps, means, sds, seed)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        f.write("# means = " + " ".join(f"{x:.6f}" for x in means) + "\n")
        f.write("# sds = " + " ".join(f"{x:.6f}" for x in sds) + "\n")
        np.savetxt(f, signal, fmt="%.6f")

    print(f"wrote {n} observations to {outfile}")
    print("true changepoints =", true_bkps)


if __name__ == "__main__":
    main()

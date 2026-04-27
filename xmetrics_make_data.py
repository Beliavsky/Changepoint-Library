"""
Create a univariate mean-shift data file for changepoint metrics comparisons.
"""

from __future__ import annotations

import sys

import numpy as np
import ruptures as rpt


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xmetrics_data.txt"
    seed = 12
    n = 1200
    n_bkps = 3
    noise_std = 1.5

    np.random.seed(seed)
    signal, true_bkps = rpt.pw_constant(
        n_samples=n,
        n_bkps=n_bkps,
        noise_std=noise_std,
    )

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write(f"# n_bkps = {n_bkps}\n")
        f.write(f"# noise_std = {noise_std}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        np.savetxt(f, signal, fmt="%.6f")

    print(f"wrote {n} observations to {outfile}")
    print("true bkps =", true_bkps)


if __name__ == "__main__":
    main()

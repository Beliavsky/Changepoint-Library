"""
Create a squared-returns data file for variance-regime RBF-PELT comparisons.
"""

import sys

import numpy as np


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xvariance_regime_data.txt"
    seed = 8
    n = 3000

    np.random.seed(seed)
    frac_bkps = np.array([0.36, 0.64])
    bkps_true = np.round(n * frac_bkps).astype(int)
    bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))
    lengths = np.diff(np.r_[0, bkps_true, n])
    sigmas = [0.005, 0.020, 0.010]

    r = np.concatenate(
        [np.random.normal(0.0, sigmas[i], lengths[i]) for i in range(len(lengths))]
    )
    x = (r**2).reshape(-1, 1)

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# {x.shape[0]} {x.shape[1]} nrow ncol\n")
        f.write(f"# seed = {seed}\n")
        f.write("# true_bkps = " + " ".join(str(v) for v in bkps_true.tolist() + [n]) + "\n")
        np.savetxt(f, x, fmt="%.6f")

    print(f"wrote {x.shape[0]} x {x.shape[1]} matrix to {outfile}")
    print("true bkps =", bkps_true.tolist() + [n])


if __name__ == "__main__":
    main()

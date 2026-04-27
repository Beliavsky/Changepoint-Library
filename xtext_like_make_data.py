"""
Create an embedding-like multivariate data file for RBF-PELT comparisons.
"""

import sys

import numpy as np


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xtext_like_data.txt"
    seed = 7
    n = 125
    p = 8

    np.random.seed(seed)
    frac_bkps = np.array([0.32, 0.60])
    bkps_true = np.round(n * frac_bkps).astype(int)
    bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))
    lengths = np.diff(np.r_[0, bkps_true, n])
    means = [0.0, 2.0, -1.0]

    emb = np.vstack(
        [np.random.normal(means[i], 1.0, size=(lengths[i], p)) for i in range(len(lengths))]
    )

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# {emb.shape[0]} {emb.shape[1]} nrow ncol\n")
        f.write(f"# seed = {seed}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in bkps_true.tolist() + [n]) + "\n")
        np.savetxt(f, emb, fmt="%.6f")

    print(f"wrote {emb.shape[0]} x {emb.shape[1]} matrix to {outfile}")
    print("true bkps =", bkps_true.tolist() + [n])


if __name__ == "__main__":
    main()

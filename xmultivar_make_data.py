"""
Create a multivariate mean-shift data file for PELT comparisons.
"""

import sys

import numpy as np


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xmultivar_data.txt"
    seed = 5
    n = 9000

    np.random.seed(seed)
    frac_bkps = np.array([1 / 3, 2 / 3])
    bkps_true = np.round(n * frac_bkps).astype(int)
    bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))
    lengths = np.diff(np.r_[0, bkps_true, n])

    means_x1 = [0.0, 2.0, -1.0]
    means_x2 = [1.0, 1.0, 4.0]
    sd_x1 = [1.0, 1.0, 1.0]
    sd_x2 = [1.0, 1.0, 1.0]

    x1 = np.concatenate(
        [np.random.normal(means_x1[i], sd_x1[i], lengths[i]) for i in range(len(lengths))]
    )
    x2 = np.concatenate(
        [np.random.normal(means_x2[i], sd_x2[i], lengths[i]) for i in range(len(lengths))]
    )
    signal = np.column_stack([x1, x2])

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# {signal.shape[0]} {signal.shape[1]} nrow ncol\n")
        f.write(f"# seed = {seed}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in bkps_true.tolist() + [n]) + "\n")
        np.savetxt(f, signal, fmt="%.6f")

    print(f"wrote {signal.shape[0]} x {signal.shape[1]} matrix to {outfile}")
    print("true bkps =", bkps_true.tolist() + [n])


if __name__ == "__main__":
    main()

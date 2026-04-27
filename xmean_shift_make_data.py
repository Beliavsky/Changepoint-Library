"""
Create a univariate mean-shift data file for PELT comparisons.
"""

import sys

import numpy as np


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xmean_shift_data.txt"
    seed = 1
    n = 10000
    means = [0.0, 3.0, -1.0]
    sds = [1.0, 1.0, 1.0]

    np.random.seed(seed)
    frac_bkps = np.array([1 / 3, 7 / 12])
    bkps_true = np.round(n * frac_bkps).astype(int)
    bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))
    lengths = np.diff(np.r_[0, bkps_true, n])

    y = np.concatenate(
        [np.random.normal(means[i], sds[i], lengths[i]) for i in range(len(lengths))]
    )

    with open(outfile, "w", encoding="ascii") as f:
        f.write(f"# seed = {seed}\n")
        f.write(f"# n = {n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in bkps_true.tolist() + [n]) + "\n")
        for val in y:
            f.write(f"{val:.6f}\n")

    print(f"wrote {n} observations to {outfile}")
    print("true bkps =", bkps_true.tolist() + [n])


if __name__ == "__main__":
    main()

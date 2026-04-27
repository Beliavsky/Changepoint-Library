#!/usr/bin/env python
"""Write shared data for changepoints::WBSIP.cov comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def write_intervals(path: Path, n: int, m: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    alpha = rng.integers(1, n + 1, size=m)
    beta = rng.integers(1, n + 1, size=m)
    lo = np.minimum(alpha, beta)
    hi = np.maximum(alpha, beta)
    with path.open("w", encoding="utf-8") as fh:
        fh.write("# random WBS intervals\n")
        for a, b in zip(lo, hi):
            fh.write(f"{int(a)} {int(b)}\n")


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xwbsip_cov_data.txt")
    prime_file = data_file.with_name(data_file.stem + "_prime.txt")
    interval_file = data_file.with_name(data_file.stem + "_intervals.txt")

    seed_x = 1234
    seed_xprime = 4321
    seed_intervals = 202
    n = 180
    p = 4
    m = 200
    true_cps = [80]

    cov1 = np.array(
        [
            [1.0, 0.65, 0.0, 0.0],
            [0.65, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.15],
            [0.0, 0.0, 0.15, 1.0],
        ]
    )
    cov2 = np.array(
        [
            [1.7, -0.15, 0.0, 0.0],
            [-0.15, 0.8, 0.0, 0.0],
            [0.0, 0.0, 0.7, -0.55],
            [0.0, 0.0, -0.55, 1.5],
        ]
    )

    rng_x = np.random.default_rng(seed_x)
    rng_xprime = np.random.default_rng(seed_xprime)

    x = np.vstack(
        (
            rng_x.multivariate_normal(np.zeros(p), cov1, size=true_cps[0]),
            rng_x.multivariate_normal(np.zeros(p), cov2, size=n - true_cps[0]),
        )
    )
    x_prime = np.vstack(
        (
            rng_xprime.multivariate_normal(np.zeros(p), cov1, size=true_cps[0]),
            rng_xprime.multivariate_normal(np.zeros(p), cov2, size=n - true_cps[0]),
        )
    )

    for path, values, seed in ((data_file, x, seed_x), (prime_file, x_prime, seed_xprime)):
        with path.open("w", encoding="utf-8") as fh:
            fh.write(f"# seed = {seed}\n")
            fh.write(f"# n = {n}\n")
            fh.write(f"# p = {p}\n")
            fh.write("# true_bkps = " + " ".join(str(v) for v in true_cps) + "\n")
            fh.write("# rows are time points and columns are variables\n")
            np.savetxt(fh, values, fmt="%.12f")

    write_intervals(interval_file, n, m, seed_intervals)

    print(f"wrote {data_file}")
    print(f"wrote {prime_file}")
    print(f"wrote {interval_file}")
    print("true changepoints =", true_cps)


if __name__ == "__main__":
    main()

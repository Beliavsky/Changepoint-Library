"""
Read a univariate series from a text file and evaluate changepoint estimates
with ruptures metrics.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import ruptures as rpt
from ruptures.metrics import hausdorff, precision_recall, randindex


def read_true_bkps(path: str) -> list[int]:
    with open(path, encoding="ascii") as f:
        for line in f:
            if line.startswith("# true_bkps ="):
                return [int(x) for x in line.split("=", 1)[1].split()]
    raise ValueError("input file is missing '# true_bkps =' metadata")


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xmetrics_data.txt"
    n_bkps = 3
    min_size = 40
    margin = 40

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)
    true_bkps = read_true_bkps(data_file)
    n = signal.shape[0]

    methods = [
        ("binseg", rpt.Binseg(model="l2", min_size=min_size, jump=1)),
        ("bottomup", rpt.BottomUp(model="l2", min_size=min_size)),
        ("window", rpt.Window(width=max(80, n // 10), model="l2")),
        ("dynp", rpt.Dynp(model="l2", min_size=min_size, jump=1)),
    ]

    print("file =", data_file)
    print("n =", n)
    print("true bkps =", true_bkps)
    print("n_bkps =", n_bkps)
    print("min_size =", min_size)
    print("margin =", margin)
    print()

    for name, algo in methods:
        t1 = time.perf_counter()
        bkps = algo.fit(signal).predict(n_bkps=n_bkps)
        elapsed = time.perf_counter() - t1
        precision, recall = precision_recall(true_bkps, bkps, margin=margin)
        hd = hausdorff(true_bkps, bkps)
        ri = randindex(true_bkps, bkps)

        print("method =", name)
        print("estimated bkps =", bkps)
        print(f"precision = {precision:.3f}")
        print(f"recall = {recall:.3f}")
        print(f"hausdorff = {hd:.3f}")
        print(f"randindex = {ri:.6f}")
        print(f"elapsed seconds = {elapsed:.3f}")
        print()

    print(f"wall time elapsed (s) = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

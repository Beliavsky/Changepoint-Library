"""
Read a series from a file and run deterministic AgglomerativeCLaP-style merging.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xagglomerative_clap import collapse_segment_process, fit_agglomerative_clap


def read_true_bkps(path: str) -> list[int]:
    """
    Read '# true_bkps = ...' metadata from the file header.
    """
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# true_bkps ="):
                return [int(tok) for tok in line.split("=", 1)[1].split()]
    raise ValueError("input file is missing '# true_bkps =' metadata")


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xagglomerative_clap_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)
    true_bkps = read_true_bkps(data_file)

    window_size = 20
    n_splits = 5
    sample_cap = 1000
    segment_labels, gain, y_true, y_pred = fit_agglomerative_clap(
        signal, true_bkps, window_size, n_splits, sample_cap
    )
    cps_out, labels_out = collapse_segment_process(true_bkps, segment_labels)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("true_bkps =", true_bkps)
    print("window_size =", window_size)
    print("n_splits =", n_splits)
    print("sample_cap =", sample_cap)
    print("segment labels =", segment_labels)
    print("collapsed bkps =", cps_out)
    print("collapsed labels =", labels_out)
    print("y_true =", y_true.tolist())
    print("y_pred =", y_pred.tolist())
    print(f"classification gain = {gain:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

"""
Read a labeled series from a file and score the CLaP-style centroid classifier.
"""

from __future__ import annotations

import sys
import time
from collections import Counter

import numpy as np

from xclap import fit_clap_centroid


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclap_data.txt"

    data = np.loadtxt(data_file, comments="#", ndmin=2)
    signal = data[:, 0].astype(float, copy=False)
    state_labels = data[:, 1].astype(int, copy=False)

    window_size = 20
    n_splits = 5
    sample_cap = 1000
    y_true, y_pred, score = fit_clap_centroid(
        signal, state_labels, window_size, n_splits, sample_cap
    )

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("window_size =", window_size)
    print("n_splits =", n_splits)
    print("sample_cap =", sample_cap)
    print("true label counts =", dict(sorted(Counter(y_true.tolist()).items())))
    print("pred label counts =", dict(sorted(Counter(y_pred.tolist()).items())))
    print("true labels =", y_true.tolist())
    print("pred labels =", y_pred.tolist())
    print(f"macro_f1 = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

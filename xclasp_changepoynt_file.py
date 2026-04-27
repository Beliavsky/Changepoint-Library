"""
Read a univariate series from a file and run changepoynt-style BinaryClaSPSegmentation.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xclasp_changepoynt import fit_clasp_changepoynt


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_changepoynt_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    detector = fit_clasp_changepoynt(signal)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("n_segments =", "learn")
    print("window_size =", "suss")
    print("n_estimators =", 10)
    print("k_neighbours =", 3)
    print("excl_radius =", 5)
    print("validation =", "significance_test")
    print("threshold =", 1e-15)
    print("estimated changepoints =", detector.predict().tolist())
    print("window_size used =", detector.window_size)
    print("n_segments used =", detector.n_segments)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

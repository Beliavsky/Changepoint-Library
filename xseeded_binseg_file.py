"""
Read a univariate series from a text file and detect changepoints with skchange.SeededBinarySegmentation.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
from skchange.change_detectors import SeededBinarySegmentation


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xseeded_binseg_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n = signal.shape[0]
    penalty = 2.0 * np.log(n)
    detector = SeededBinarySegmentation(
        penalty=penalty,
        max_interval_length=max(200, n // 8),
        growth_factor=1.5,
        selection_method="greedy",
    )
    result = detector.fit_predict(pd.DataFrame({"x": signal}))

    print("file =", data_file)
    print("n =", n)
    print(f"penalty = {penalty:.6f}")
    print("estimated changepoints =", result["ilocs"].to_list())
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

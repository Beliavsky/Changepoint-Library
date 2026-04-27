"""
Read a univariate series from a text file and detect changepoints with skchange.CROPS.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
from skchange.change_detectors import CROPS
from skchange.costs import L2Cost


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcrops_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    detector = CROPS(
        cost=L2Cost(),
        min_penalty=5.0,
        max_penalty=60.0,
        selection_method="bic",
        min_segment_length=10,
        step_size=1,
        split_cost=0.0,
        prune=True,
        pruning_margin=0.0,
    )
    result = detector.fit_predict(pd.DataFrame({"x": signal}))

    print("file =", data_file)
    print("n =", signal.shape[0])
    print(f"optimal penalty = {float(detector.optimal_penalty):.6f}")
    print("path solutions =", len(detector.change_points_lookup))
    print("estimated changepoints =", result["ilocs"].to_list())
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

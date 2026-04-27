"""
Read a series from a file and score a deterministic roerich RuLSIF-style detector.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xroerich_rulsif import fit_roerich_rulsif


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xroerich_rulsif_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    window_size = 40
    periods = 1
    step = 1
    alpha = 0.05
    l2 = 1e-3
    cp, score = fit_roerich_rulsif(signal, window_size, periods, step, 2357, alpha, l2)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("window_size =", window_size)
    print("periods =", periods)
    print("step =", step)
    print("base_regressor = linear_rulsif")
    print("metric = pesym")
    print("alpha =", alpha)
    print("l2 =", l2)
    print("estimated changepoint =", cp)
    print(f"max raw score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

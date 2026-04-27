"""
Read a data file and score it with a deterministic Gaussian-kernel ULSIF detector.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from xrulsif_utils import fit_rulsif_gaussian
from xsst_utils import read_series_file


def main() -> None:
    t0 = time.perf_counter()
    data_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xulsif_data.txt")
    signal, true_cp = read_series_file(data_path)

    window_length = 10
    n_windows = 50
    lag = 50
    sigma = 2.0
    lambda_reg = 1e-2
    scoring_step = 1

    _, cp, best_score = fit_rulsif_gaussian(
        signal,
        window_length=window_length,
        n_windows=n_windows,
        lag=lag,
        alpha=0.0,
        sigma=sigma,
        lambda_reg=lambda_reg,
        scoring_step=scoring_step,
        symmetric=True,
    )

    print("file =", data_path)
    print("n =", signal.size)
    if true_cp is not None:
        print("true changepoint =", true_cp)
    print("window_length =", window_length)
    print("n_windows =", n_windows)
    print("lag =", lag)
    print("alpha = 0.0")
    print("sigma =", sigma)
    print("lambda =", lambda_reg)
    print("scoring_step =", scoring_step)
    print("estimated changepoint =", cp)
    print(f"max raw score = {best_score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

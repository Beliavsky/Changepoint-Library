"""
Read a data file and score it with a deterministic ESST-style detector.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

from xesst_utils import fit_esst_exact, read_series_file


def main() -> None:
    t0 = time.perf_counter()
    data_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xesst_data.txt")
    signal, true_cp = read_series_file(data_path)

    window_length = 48
    n_windows = 24
    lag = 24
    rank = 2
    scoring_step = 1

    _, cp, best_score = fit_esst_exact(
        signal,
        window_length=window_length,
        n_windows=n_windows,
        lag=lag,
        rank=rank,
        scoring_step=scoring_step,
    )

    print("file =", data_path)
    print("n =", signal.size)
    if true_cp is not None:
        print("true changepoint =", true_cp)
    print("window_length =", window_length)
    print("n_windows =", n_windows)
    print("lag =", lag)
    print("rank =", rank)
    print("scoring_step =", scoring_step)
    print("method = exact")
    print("estimated changepoint =", cp)
    print(f"max raw score = {best_score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

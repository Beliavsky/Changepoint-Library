"""
Read a univariate series from a text file and score it with changepoynt.BOCPD.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xbocpd_changepoynt_utils import fit_changepoynt_bocpd


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xbocpd_changepoynt_data.txt"
    run_length = 120

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    _, cp_raw, cp_post, best_score, prior_mean, prior_var, signal_var = fit_changepoynt_bocpd(signal, run_length)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("run_length =", run_length)
    print("raw score argmax =", cp_raw)
    print("post_warmup argmax =", cp_post)
    print(f"max score = {best_score:.6f}")
    print(f"prior_mean = {prior_mean:.6f}")
    print(f"prior_var = {prior_var:.6f}")
    print(f"signal_var = {signal_var:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

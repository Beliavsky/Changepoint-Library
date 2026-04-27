"""
Simulate a univariate series and score it with changepoynt.BOCPD.
"""

from __future__ import annotations

import time

from xbocpd_changepoynt_utils import fit_changepoynt_bocpd, simulate_bocpd_signal


def main() -> None:
    t0 = time.perf_counter()
    run_length = 120
    true_bkps = [0, 200, 400]
    means = [0.0, 3.5, -1.5]
    sds = [1.0, 1.0, 1.6]
    signal = simulate_bocpd_signal(600, true_bkps, means, sds, 21)

    _, cp_raw, cp_post, best_score, prior_mean, prior_var, signal_var = fit_changepoynt_bocpd(signal, run_length)

    print("n =", signal.shape[0])
    print("run_length =", run_length)
    print("true changepoints =", true_bkps)
    print("raw score argmax =", cp_raw)
    print("post_warmup argmax =", cp_post)
    print(f"max score = {best_score:.6f}")
    print(f"prior_mean = {prior_mean:.6f}")
    print(f"prior_var = {prior_var:.6f}")
    print(f"signal_var = {signal_var:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

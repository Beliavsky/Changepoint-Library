"""
Shared helpers for changepoynt-specific BOCPD comparison scripts.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CHANGEPOYNT_REPO = Path(r"c:/python/public_domain/github/changepoynt")
if str(CHANGEPOYNT_REPO) not in sys.path:
    sys.path.insert(0, str(CHANGEPOYNT_REPO))

from changepoynt.algorithms.bocpd import BOCPD


def simulate_bocpd_signal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with piecewise mean and variance.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    for i, start in enumerate(regime_starts):
        end = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        signal[start:end] = rng.normal(loc=means[i], scale=sds[i], size=end - start)
    return signal


def write_series_file(path: str | Path, signal: np.ndarray, true_bkps: list[int], run_length: int) -> None:
    """
    Write a one-column text file with BOCPD metadata in header comments.
    """
    path = Path(path)
    with path.open("w", encoding="ascii") as handle:
        handle.write(f"# run_length = {run_length}\n")
        handle.write("# true_bkps = " + " ".join(str(x) for x in true_bkps) + "\n")
        np.savetxt(handle, signal, fmt="%.6f")


def fit_changepoynt_bocpd(
    signal: np.ndarray,
    run_length: int,
    change_length_threshold: int | None = None,
) -> tuple[np.ndarray, int, int, float, float, float, float]:
    """
    Run changepoynt BOCPD and return score plus raw and post-warmup argmax summaries.
    """
    detector = BOCPD(run_length=run_length, change_length_threshold=change_length_threshold)
    score = np.asarray(detector.transform(signal), dtype=float)
    cp_raw = int(np.argmax(score))
    cp_post = int(np.argmax(score[run_length - 1 :]) + run_length - 1)
    return (
        score,
        cp_raw,
        cp_post,
        float(score[cp_raw]),
        float(detector.prior_mean),
        float(detector.prior_var),
        float(detector.signal_var),
    )

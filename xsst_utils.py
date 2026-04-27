"""
Shared helpers for changepoynt SST comparison scripts.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CHANGEPOYNT_REPO = Path(r"c:/python/public_domain/github/changepoynt")
if str(CHANGEPOYNT_REPO) not in sys.path:
    sys.path.insert(0, str(CHANGEPOYNT_REPO))

from changepoynt.algorithms.sst import SST


def simulate_frequency_change_signal(
    n_per_segment: int,
    period_before: int,
    period_after: int,
    noise: float,
    seed: int,
) -> tuple[np.ndarray, int]:
    """
    Simulate the frequency-change signal used for the SST compare case.
    """
    rng = np.random.default_rng(seed)
    t = np.arange(n_per_segment)
    left = np.sin(2.0 * np.pi * t / period_before)
    right = np.sin(2.0 * np.pi * t / period_after)
    signal = np.concatenate([left, right]).astype(float)
    signal += noise * rng.standard_normal(signal.shape[0])
    return signal, n_per_segment


def write_series_file(path: str | Path, signal: np.ndarray, true_cp: int) -> None:
    """
    Write a one-column text file with a header comment for the true changepoint.
    """
    path = Path(path)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"# true changepoint = {true_cp}\n")
        for value in signal:
            handle.write(f"{value:.17g}\n")


def read_series_file(path: str | Path) -> tuple[np.ndarray, int | None]:
    """
    Read the one-column text format and return the optional true changepoint from the header.
    """
    path = Path(path)
    true_cp: int | None = None
    values: list[float] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            lower = line.lower()
            if "true changepoint" in lower and "=" in line:
                try:
                    true_cp = int(line.split("=", 1)[1].strip())
                except ValueError:
                    true_cp = None
            continue
        values.append(float(line))
    return np.asarray(values, dtype=float), true_cp


def fit_sst_naive(
    signal: np.ndarray,
    window_length: int,
    n_windows: int,
    lag: int,
    rank: int,
    scoring_step: int = 1,
    scale: bool = True,
) -> tuple[np.ndarray, int, float]:
    """
    Run changepoynt SST with the exact slow naive method and return the full score vector and argmax.
    """
    detector = SST(
        window_length=window_length,
        n_windows=n_windows,
        lag=lag,
        rank=rank,
        scale=scale,
        method="naive",
        scoring_step=scoring_step,
        use_fast_hankel=False,
    )
    score = np.asarray(detector.transform(signal), dtype=float)
    cp = int(np.argmax(score))
    return score, cp, float(score[cp])

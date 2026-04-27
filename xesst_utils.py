"""
Shared helpers for deterministic ESST-style comparison scripts.
"""

from __future__ import annotations

import numpy as np

from xsst_utils import read_series_file, simulate_frequency_change_signal, write_series_file


def fit_esst_exact(
    signal: np.ndarray,
    window_length: int,
    n_windows: int,
    lag: int,
    rank: int,
    scoring_step: int = 1,
    scale: bool = True,
) -> tuple[np.ndarray, int, float]:
    """
    Run a deterministic exact-SVD ESST analog using the package's score definition.
    """
    x = np.asarray(signal, dtype=float)
    if scale:
        xmin = float(np.min(x))
        xmax = float(np.max(x))
        if xmax == xmin:
            x = x - xmin
        else:
            x = (x - xmin) / (xmax - xmin)
        x = x + 1.0
    else:
        x = x.copy()

    score = np.zeros_like(x)
    start_idx = window_length + n_windows + lag
    offset = n_windows + lag

    for idx in range(start_idx, x.shape[0], scoring_step):
        hankel_past = np.empty((window_length, n_windows), dtype=float)
        hankel_future = np.empty((window_length, n_windows), dtype=float)
        for cx in range(n_windows):
            col = n_windows - cx - 1
            hankel_past[:, col] = x[idx - lag - window_length - cx : idx - lag - cx]
            hankel_future[:, col] = x[idx - window_length - cx : idx - cx]
        hankel = np.concatenate((hankel_past, hankel_future), axis=1)

        _, singular_values, vt = np.linalg.svd(hankel, full_matrices=False)
        vt = vt[:rank, :]
        singular_values = singular_values[:rank]
        vt = vt - np.min(vt, axis=1, keepdims=True) + 1.0
        vt = vt / np.sum(vt, axis=1, keepdims=True)
        half = vt.shape[1] // 2
        skew = np.abs(np.mean(vt[:, :half] - vt[:, half:], axis=1))
        local_score = float(np.dot(singular_values, skew) / np.sum(singular_values))

        left_idx = max(idx - offset - scoring_step // 2, 0)
        right_idx = min(idx - offset + (scoring_step + 1) // 2, x.shape[0])
        score[left_idx:right_idx] = local_score

    cp = int(np.argmax(score))
    return score, cp, float(score[cp])

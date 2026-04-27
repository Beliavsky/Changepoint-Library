"""
Shared Python helpers for change-point comparison scripts.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


def build_prefix_1d(signal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Build prefix sums for a univariate signal.
    """
    x = np.ravel(np.asarray(signal, dtype=float))
    sx = np.empty(x.size + 1, dtype=float)
    sxx = np.empty(x.size + 1, dtype=float)
    sx[0] = 0.0
    sxx[0] = 0.0
    np.cumsum(x, out=sx[1:])
    np.cumsum(x * x, out=sxx[1:])
    return sx, sxx


def segment_sse_1d(
    sx: np.ndarray,
    sxx: np.ndarray,
    start: int,
    end: int,
    min_size: int,
) -> float:
    """
    Compute univariate l2 segment SSE on the half-open interval [start, end).
    """
    length = end - start
    if length < min_size:
        return float("inf")
    sum_x = sx[end] - sx[start]
    sum_xx = sxx[end] - sxx[start]
    return float(sum_xx - (sum_x * sum_x) / length)


def refine_bkps_l2_1d(
    signal: np.ndarray,
    bkps: Iterable[int],
    min_size: int,
    max_iter: int = 20,
) -> list[int]:
    """
    Refine internal breakpoints by local l2 re-optimization with neighbors fixed.
    """
    refined = [int(b) for b in bkps]
    if not refined:
        return refined

    n = np.ravel(np.asarray(signal)).size
    if refined[-1] != n:
        raise ValueError("breakpoints must include the final endpoint n")

    sx, sxx = build_prefix_1d(signal)

    for _ in range(max_iter):
        changed = False
        for i in range(len(refined) - 1):
            left = 0 if i == 0 else refined[i - 1]
            right = refined[i + 1]
            lo = left + min_size
            hi = right - min_size
            if lo > hi:
                continue

            best_k = refined[i]
            best_cost = (
                segment_sse_1d(sx, sxx, left, best_k, min_size)
                + segment_sse_1d(sx, sxx, best_k, right, min_size)
            )

            for k in range(lo, hi + 1):
                cost = (
                    segment_sse_1d(sx, sxx, left, k, min_size)
                    + segment_sse_1d(sx, sxx, k, right, min_size)
                )
                if cost < best_cost:
                    best_cost = cost
                    best_k = k

            if best_k != refined[i]:
                refined[i] = best_k
                changed = True

        if not changed:
            break

    return refined

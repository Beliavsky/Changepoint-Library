"""
Simulate a univariate series and score a deterministic ClaSPEnsemble-style split.

This mirrors claspy.clasp.ClaSPEnsemble, but uses a shared deterministic temporal
constraint generator so Python and Fortran can be compared exactly.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

CLASPY_REPO = Path(r"c:/python/public_domain/github/claspy")
if str(CLASPY_REPO) not in sys.path:
    sys.path.insert(0, str(CLASPY_REPO))

from claspy.clasp import ClaSP

from xclasp import simulate_clasp_signal


def park_miller_next(state: int) -> int:
    """
    Advance the Park-Miller generator by one step.
    """
    a = 16807
    m = 2147483647
    q = 127773
    r = 2836
    hi = state // q
    lo = state % q
    test = a * lo - r * hi
    return test if test > 0 else test + m


def temporal_constraints(
    n: int,
    n_estimators: int,
    window_size: int,
    excl_radius: int,
    seed: int,
) -> list[tuple[int, int]]:
    """
    Generate deterministic temporal constraints in the same style as ClaSPEnsemble.
    """
    min_seg_size = window_size * excl_radius
    constraints: list[tuple[int, int]] = [(0, n)]
    state = abs(seed) % 2147483646 + 1

    while len(constraints) < n_estimators and n > 3 * min_seg_size:
        state = park_miller_next(state)
        lbound = (state - 1) % n
        state = park_miller_next(state)
        area = (state - 1) % n
        if n - lbound < area:
            area = n - lbound
        ubound = lbound + area
        if ubound - lbound < 2 * min_seg_size:
            continue
        constraints.append((lbound, ubound))

    constraints.sort(key=lambda tc: tc[1] - tc[0], reverse=True)
    return constraints


def fit_clasp_ensemble(
    signal: np.ndarray,
    n_estimators: int,
    window_size: int,
    k_neighbours: int,
    excl_radius: int,
    score: str = "roc_auc",
    early_stopping: bool = True,
    seed: int = 2357,
) -> tuple[int | None, float, list[tuple[int, int]]]:
    """
    Run deterministic ClaSPEnsemble-style interval selection and return the best split.
    """
    constraints = temporal_constraints(
        signal.shape[0], n_estimators, window_size, excl_radius, seed
    )
    best_cp: int | None = None
    best_score = -np.inf

    for idx, (lbound, ubound) in enumerate(constraints):
        clasp = ClaSP(
            window_size=window_size,
            k_neighbours=k_neighbours,
            distance="znormed_euclidean_distance",
            score=score,
            excl_radius=excl_radius,
            n_jobs=1,
        )
        profile = clasp.fit_transform(signal[lbound:ubound])
        profile = (profile + (ubound - lbound) / signal.shape[0]) / 2.0
        local_score = float(np.nanmax(profile))
        local_cp = lbound + int(np.nanargmax(profile))

        if local_score > best_score or (best_cp is None and idx == len(constraints) - 1):
            best_score = local_score
            best_cp = local_cp
        else:
            if early_stopping:
                break

    return best_cp, best_score, constraints


def main() -> None:
    t0 = time.perf_counter()
    seed = 101
    n = 240
    true_cp = 120
    n_estimators = 5
    window_size = 12
    k_neighbours = 3
    excl_radius = 5

    signal = simulate_clasp_signal(n, [0, true_cp], [0.0, 4.0], [1.0, 1.0], seed)
    cp, score, constraints = fit_clasp_ensemble(
        signal, n_estimators, window_size, k_neighbours, excl_radius
    )

    print("n =", n)
    print("true changepoint =", true_cp)
    print("n_estimators =", n_estimators)
    print("window_size =", window_size)
    print("k_neighbours =", k_neighbours)
    print("excl_radius =", excl_radius)
    print("temporal constraints =", constraints)
    print("estimated changepoint =", cp)
    print(f"max profile score = {score:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

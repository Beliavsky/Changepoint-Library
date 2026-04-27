"""
Read a univariate series from a file and run the deterministic ClaSPEnsemble analog.
"""

from __future__ import annotations

import sys
import time

import numpy as np

from xclasp_ensemble import fit_clasp_ensemble


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclasp_ensemble_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    n_estimators = 5
    window_size = 12
    k_neighbours = 3
    excl_radius = 5
    cp, score, constraints = fit_clasp_ensemble(
        signal, n_estimators, window_size, k_neighbours, excl_radius
    )

    print("file =", data_file)
    print("n =", signal.shape[0])
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

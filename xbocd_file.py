#!/usr/bin/env python
"""Run the local `bocd` package on a one-column data file."""

from __future__ import annotations

import sys
import time

from xbocd_utils import read_series_file, run_bocd_student_t


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xbocd_data.txt"
    lam = 100.0

    x = read_series_file(data_file)
    t0 = time.perf_counter()
    out = run_bocd_student_t(x, lam=lam, mu=0.0, kappa=1.0, alpha=1.0, beta=1.0)
    elapsed = time.perf_counter() - t0

    print(f"file = {data_file}")
    print(f"n = {x.size}")
    print(f"lambda = {lam}")
    print("distribution = StudentT")
    print("hazard = ConstantHazard")
    print(f"estimated changepoints = {out['map_cps']}")
    print(f"final MAP run length = {out['final_rt']}")
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

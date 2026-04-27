#!/usr/bin/env python
"""Fit a regular multi-series BEAST123-style model from a matrix data file."""

from __future__ import annotations

import sys
import time

from xbeast_utils import read_matrix_file, solve_beast123_regular, top_cp_probabilities


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xbeast123_data.txt"
    max_cp = 4
    min_seg_len = 30

    x = read_matrix_file(data_file)
    t0 = time.perf_counter()
    out = solve_beast123_regular(x, max_cp=max_cp, min_seg_len=min_seg_len)
    elapsed = time.perf_counter() - t0

    print(f"file = {data_file}")
    print(f"n = {x.shape[0]}")
    print(f"p = {x.shape[1]}")
    print("season = none")
    print("beast123 regular interface")
    print(f"max_cp = {max_cp}")
    print(f"min_seg_len = {min_seg_len}")
    for j, cps in enumerate(out["best_cps"], start=1):
        print(f"series {j} changepoints = {cps}")
    print("top mean cp probabilities =", [(i, round(pv, 6)) for i, pv in top_cp_probabilities(out["mean_cp_prob"])])
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Fit a regular trend-only BEAST-style model from a one-column data file."""

from __future__ import annotations

import sys
import time

from xbeast_utils import read_series_file, solve_beast_trend_only, top_cp_probabilities


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xbeast_data.txt"
    max_cp = 4
    min_seg_len = 30

    y = read_series_file(data_file)
    t0 = time.perf_counter()
    out = solve_beast_trend_only(y, max_cp=max_cp, min_seg_len=min_seg_len)
    elapsed = time.perf_counter() - t0

    print(f"file = {data_file}")
    print(f"n = {y.size}")
    print("season = none")
    print(f"max_cp = {max_cp}")
    print(f"min_seg_len = {min_seg_len}")
    print(f"estimated changepoints = {out['best_cps']}")
    print("model weights =", [round(float(w), 6) for w in out["model_weights"]])
    print("top cp probabilities =", [(i, round(p, 6)) for i, p in top_cp_probabilities(out["cp_prob"])])
    print(f"elapsed seconds = {elapsed:.3f}")


if __name__ == "__main__":
    main()

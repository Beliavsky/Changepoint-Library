#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xbcp_multivar_data.txt")
    cps = [40, 80]
    rng = np.random.default_rng(123)
    y = np.column_stack(
        [
            np.concatenate(
                [
                    rng.normal(0.0, 0.6, size=40),
                    rng.normal(2.5, 0.6, size=40),
                    rng.normal(-1.2, 0.6, size=40),
                ]
            ),
            np.concatenate(
                [
                    rng.normal(1.0, 0.7, size=40),
                    rng.normal(-1.5, 0.7, size=40),
                    rng.normal(2.2, 0.7, size=40),
                ]
            ),
        ]
    )
    np.savetxt(out_path, y, fmt="%.12f", header=f"true_bkps = {' '.join(map(str, cps))}")
    print(f"wrote {out_path.name} with shape {y.shape[0]} x {y.shape[1]}")
    print(f"true changepoints = {cps}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xsegmented_multi_data.txt")
    rng = np.random.default_rng(456)

    x = np.arange(1, 181, dtype=float)
    y = np.where(
        x < 60.0,
        2.0 + 0.15 * x,
        np.where(
            x < 120.0,
            2.0 + 0.15 * 60.0 + 0.55 * (x - 60.0),
            2.0 + 0.15 * 60.0 + 0.55 * 60.0 - 0.10 * (x - 120.0),
        ),
    ) + rng.normal(0.0, 1.2, size=x.size)

    data = np.column_stack([x, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {x.size} rows")
    print("true breakpoints = [60, 120]")


if __name__ == "__main__":
    main()

#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xsegmented_data.txt")
    rng = np.random.default_rng(123)

    x = np.arange(1, 121, dtype=float)
    y = np.where(
        x < 60.0,
        5.0 + 0.2 * x,
        5.0 + 0.2 * 60.0 + 0.8 * (x - 60.0),
    ) + rng.normal(0.0, 1.5, size=x.size)

    data = np.column_stack([x, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {x.size} rows")
    print("true breakpoint = 60")


if __name__ == "__main__":
    main()

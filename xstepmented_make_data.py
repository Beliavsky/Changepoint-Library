#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xstepmented_data.txt")
    rng = np.random.default_rng(321)

    x = np.arange(1, 121, dtype=float)
    y = np.where(x < 60.0, 4.0, 12.0) + rng.normal(0.0, 1.0, size=x.size)

    data = np.column_stack([x, y])
    np.savetxt(out, data, fmt="%.10f")

    print(f"wrote {out} with {x.size} rows")
    print("true jumpoint = 60")


if __name__ == "__main__":
    main()

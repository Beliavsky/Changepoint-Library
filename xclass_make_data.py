"""
Write the simulated ClaSS example to a text file.
"""

from __future__ import annotations

import sys
from pathlib import Path

from xclass import simulate_class_signal


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xclass_data.txt")

    seed = 101
    n = 240
    true_cp = 120
    means = (0.0, 4.0)
    sds = (1.0, 1.0)
    signal = simulate_class_signal(n, true_cp, means, sds, seed)

    with data_file.open("w", encoding="utf-8") as fh:
        fh.write(f"# n = {n}\n")
        fh.write(f"# true changepoint = {true_cp}\n")
        fh.write(f"# seed = {seed}\n")
        for value in signal:
            fh.write(f"{value:.16f}\n")

    print(f"wrote {n} observations to {data_file}")
    print(f"true changepoint = {true_cp}")


if __name__ == "__main__":
    main()

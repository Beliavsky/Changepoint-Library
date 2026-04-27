"""
Create a data file for the deterministic RuLSIF comparison workflow.
"""

from __future__ import annotations

import sys
from pathlib import Path

from xrulsif_utils import simulate_rulsif_signal
from xsst_utils import write_series_file


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xrulsif_data.txt")
    true_cp = 200
    signal = simulate_rulsif_signal(400, true_cp, (0.0, 3.0), (1.0, 1.0), 101)
    write_series_file(out_path, signal, true_cp)
    print(f"wrote {signal.size} observations to {out_path}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

"""
Create a data file for the deterministic ESST comparison workflow.
"""

from __future__ import annotations

import sys
from pathlib import Path

from xesst_utils import simulate_frequency_change_signal, write_series_file


def main() -> None:
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xesst_data.txt")
    n_per_segment = 320
    true_cp = n_per_segment
    signal, _ = simulate_frequency_change_signal(
        n_per_segment=n_per_segment,
        period_before=48,
        period_after=14,
        noise=0.02,
        seed=1234,
    )
    write_series_file(out_path, signal, true_cp)
    print(f"wrote {signal.size} observations to {out_path}")
    print("true changepoint =", true_cp)


if __name__ == "__main__":
    main()

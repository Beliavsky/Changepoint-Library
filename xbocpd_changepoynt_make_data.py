"""
Create a univariate Gaussian data file for changepoynt.BOCPD comparisons.
"""

from __future__ import annotations

import sys

from xbocpd_changepoynt_utils import simulate_bocpd_signal, write_series_file


def main() -> None:
    outfile = sys.argv[1] if len(sys.argv) > 1 else "xbocpd_changepoynt_data.txt"
    run_length = 120
    true_bkps = [0, 200, 400]
    means = [0.0, 3.5, -1.5]
    sds = [1.0, 1.0, 1.6]
    signal = simulate_bocpd_signal(600, true_bkps, means, sds, 21)
    write_series_file(outfile, signal, true_bkps, run_length)
    print(f"wrote {signal.shape[0]} observations to {outfile}")
    print("true changepoints =", true_bkps)


if __name__ == "__main__":
    main()

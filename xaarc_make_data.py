"""
Write contaminated input data for changepoints::aARC comparisons.
"""

from __future__ import annotations

import sys

from xarc_make_data import simulate_arc_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xaarc_data.txt"
    y, cps, outlier_idx, outlier_shift = simulate_arc_signal(seed=0)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# contaminated univariate mean-change data for changepoints::aARC\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        fh.write("# outlier_idx = " + " ".join(str(i) for i in outlier_idx) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {len(y)} observations to {data_file}")
    print("true changepoints =", cps)
    print("outlier_idx =", outlier_idx)
    print("outlier_shift =", outlier_shift)


if __name__ == "__main__":
    main()

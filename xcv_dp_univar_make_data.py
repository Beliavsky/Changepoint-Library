"""
Write simulated input data for changepoints::CV.search.DP.univar comparisons.
"""

from __future__ import annotations

import sys

from xdp_univar_make_data import simulate_dp_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xcv_dp_univar_data.txt"
    y, cps = simulate_dp_signal(seed=0)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# univariate mean-change data for changepoints::CV.search.DP.univar\n")
        fh.write("# true_cps = " + " ".join(str(cp) for cp in cps) + "\n")
        for value in y:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {len(y)} observations to {data_file}")
    print("true changepoints =", cps)


if __name__ == "__main__":
    main()

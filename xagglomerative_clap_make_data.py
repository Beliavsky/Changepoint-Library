"""
Write repeated-state input data for AgglomerativeCLaP-style comparisons.
"""

from __future__ import annotations

import sys

from xclasp import simulate_clasp_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xagglomerative_clap_data.txt"

    n = 480
    true_bkps = [80, 160, 240, 320, 400]
    true_states = [1, 2, 3, 1, 2, 3]
    means = [0.0, 2.0, -2.0, 0.0, 2.0, -2.0]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    signal = simulate_clasp_signal(n, [0] + true_bkps, means, sds, 101)

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# repeated-state data for AgglomerativeCLaP-style detection\n")
        fh.write("# true_bkps = " + " ".join(map(str, true_bkps)) + "\n")
        fh.write("# true_states = " + " ".join(map(str, true_states)) + "\n")
        for value in signal:
            fh.write(f"{value:.12f}\n")

    print(f"wrote {n} observations to {data_file}")
    print("true_bkps =", true_bkps)
    print("true states =", true_states)


if __name__ == "__main__":
    main()

"""
Write simulated input data for the CLaP-style state-classification comparison.
"""

from __future__ import annotations

import sys

from xclap import simulate_clap_signal


def main() -> None:
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xclap_data.txt"

    seed = 211
    n = 480
    regime_starts = [0, 80, 160, 240, 320, 400]
    phis = [0.8, -0.6, 0.2, 0.8, -0.6, 0.2]
    sigmas = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
    regime_labels = [1, 2, 3, 1, 2, 3]
    signal, state_labels = simulate_clap_signal(
        n, regime_starts, phis, sigmas, regime_labels, seed
    )

    with open(data_file, "w", encoding="utf-8") as fh:
        fh.write("# two-column AR(1) data for a CLaP-style state-classification comparison\n")
        fh.write("# regime_starts = " + " ".join(map(str, regime_starts)) + "\n")
        fh.write("# regime_labels = " + " ".join(map(str, regime_labels)) + "\n")
        fh.write("# phis = " + " ".join(map(str, phis)) + "\n")
        for value, label in zip(signal, state_labels, strict=True):
            fh.write(f"{value:.12f} {label:d}\n")

    print(f"wrote {n} labeled observations to {data_file}")
    print("regime starts =", regime_starts)
    print("regime labels =", regime_labels)
    print("phis =", phis)


if __name__ == "__main__":
    main()

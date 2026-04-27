"""
Read a count series from a text file and detect changepoints with changepoint.Bocpd.
"""

from __future__ import annotations

import sys
import time

import changepoint as cpt
import numpy as np


def main() -> None:
    t0 = time.perf_counter()
    data_file = sys.argv[1] if len(sys.argv) > 1 else "xbocpd_poisson_data.txt"

    signal = np.loadtxt(data_file, comments="#", ndmin=1, dtype=np.int64)
    if signal.ndim != 1:
        signal = np.ravel(signal)

    bocpd = cpt.Bocpd(cpt.PoissonGamma(), 120.0)
    rs = [bocpd.step(int(x)) for x in signal]
    map_cps = cpt.map_changepoints(rs)

    print("file =", data_file)
    print("n =", signal.shape[0])
    print("estimated changepoints =", map_cps)
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

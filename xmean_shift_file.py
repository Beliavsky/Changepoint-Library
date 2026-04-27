"""
Read a univariate series from a text file and detect mean shifts using PELT with the l2 cost.
"""

import sys
import time

import numpy as np
import ruptures as rpt

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xmean_shift_data.txt"
pen = 10.0

signal = np.loadtxt(data_file, comments="#", ndmin=1)
if signal.ndim != 1:
    signal = np.ravel(signal)

n = signal.shape[0]
min_size = max(5, n // 20)

algo = rpt.Pelt(model="l2", min_size=min_size).fit(signal)
bkps_est = algo.predict(pen=pen)

print("file =", data_file)
print("n =", n)
print("pen =", pen)
print("min_size =", min_size)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

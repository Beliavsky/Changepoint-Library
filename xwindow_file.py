"""
Read a univariate series from a text file and detect mean shifts using the sliding-window l2 method.
"""

import sys
import time

import numpy as np
import ruptures as rpt

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xwindow_data.txt"
n_bkps = 3

signal = np.loadtxt(data_file, comments="#", ndmin=1)
if signal.ndim != 1:
    signal = np.ravel(signal)

n = signal.shape[0]
width = max(20, n // 10)

algo = rpt.Window(width=width, model="l2").fit(signal)
bkps_est = algo.predict(n_bkps=n_bkps)

print("file =", data_file)
print("n =", n)
print("width =", width)
print("n_bkps =", n_bkps)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

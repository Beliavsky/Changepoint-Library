"""
Read a univariate series from a text file and detect a fixed number of mean shifts.
Algorithm: bottom-up segmentation with the l2 cost.
"""

import sys
import time

import numpy as np
import ruptures as rpt
from utils_py import refine_bkps_l2_1d

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xbottomup_data.txt"
n_bkps = 4

signal = np.loadtxt(data_file, comments="#", ndmin=1)
if signal.ndim != 1:
    signal = np.ravel(signal)

n = signal.shape[0]
min_size = max(5, n // 25)

algo = rpt.BottomUp(model="l2", min_size=min_size, jump=1).fit(signal)
bkps_init = algo.predict(n_bkps=n_bkps)
bkps_est = refine_bkps_l2_1d(signal, bkps_init, min_size)

print("file =", data_file)
print("n =", n)
print("n_bkps =", n_bkps)
print("min_size =", min_size)
print("initial bkps =", bkps_init)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

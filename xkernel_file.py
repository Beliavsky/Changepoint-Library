"""
Read a univariate series from a text file and detect changepoints
using kernel dynamic programming with an RBF kernel.
"""

import sys
import time

import numpy as np
import ruptures as rpt

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xkernel_data.txt"
n_bkps = 3

signal = np.loadtxt(data_file, comments="#", ndmin=1)
if signal.ndim != 1:
    signal = np.ravel(signal)

n = signal.shape[0]

algo = rpt.KernelCPD(kernel="rbf").fit(signal)
bkps_est = algo.predict(n_bkps=n_bkps)

print("file =", data_file)
print("n =", n)
print("n_bkps =", n_bkps)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

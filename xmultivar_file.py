"""
Read a multivariate series from a text file and detect common mean-shift changepoints using PELT.
"""

import sys
import time

import numpy as np
import ruptures as rpt

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xmultivar_data.txt"
pen = 12.0

signal = np.loadtxt(data_file, comments="#", ndmin=2)
n = signal.shape[0]
min_size = max(5, n // 20)

algo = rpt.Pelt(model="l2", min_size=min_size).fit(signal)
bkps_est = algo.predict(pen=pen)

print("file =", data_file)
print("n =", n)
print("p =", signal.shape[1])
print("pen =", pen)
print("min_size =", min_size)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

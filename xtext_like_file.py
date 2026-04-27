"""
Read an embedding-like matrix from a text file and detect distribution shifts using PELT with RBF cost.
"""

import sys
import time

import numpy as np
import ruptures as rpt

t0 = time.perf_counter()

data_file = sys.argv[1] if len(sys.argv) > 1 else "xtext_like_data.txt"
pen = 15.0

emb = np.loadtxt(data_file, comments="#", ndmin=2)
n, p = emb.shape
min_size = max(5, n // 20)

algo = rpt.Pelt(model="rbf", min_size=min_size).fit(emb)
bkps_est = algo.predict(pen=pen)

print("file =", data_file)
print("n =", n)
print("p =", p)
print("pen =", pen)
print("min_size =", min_size)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

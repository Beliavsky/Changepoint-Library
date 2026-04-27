"""
Detect a fixed number of mean shifts in a univariate series.
Algorithm: bottom-up segmentation with the l2 cost.
"""

import time

import numpy as np
import ruptures as rpt
from utils_py import refine_bkps_l2_1d

t0 = time.perf_counter()
np.random.seed(8)

write_data = False
outfile = "xbottomup_data.txt"

n = 100000
n_bkps = 4
noise_std = 2.0

signal, true_bkps = rpt.pw_constant(
    n_samples=n,
    n_bkps=n_bkps,
    noise_std=noise_std,
)

if write_data:
    np.savetxt(outfile, signal, fmt="%.6f")

min_size = max(5, n // 25)
algo = rpt.BottomUp(model="l2", min_size=min_size).fit(signal)
bkps_init = algo.predict(n_bkps=n_bkps)
bkps_est = refine_bkps_l2_1d(signal, bkps_init, min_size)

print("n =", n)
print("true bkps =", true_bkps)
print("initial bkps =", bkps_init)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

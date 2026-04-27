"""
Detect mean shifts in a univariate series using local window comparisons.
Algorithm: sliding-window method with the l2 cost.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(4)

write_data = False
outfile = "xwindow_data.txt"

n = 400000
n_bkps = 3
noise_std = 1.5

signal, true_bkps = rpt.pw_constant(
    n_samples=n,
    n_bkps=n_bkps,
    noise_std=noise_std
)

if write_data:
    np.savetxt(outfile, signal, fmt="%.6f")

width = max(20, n // 10)
algo = rpt.Window(width=width, model="l2").fit(signal)
bkps_est = algo.predict(n_bkps=n_bkps)

print("n =", n)
print("width =", width)
print("true bkps =", true_bkps)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

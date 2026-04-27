"""
Detect a fixed number of mean shifts in a univariate series.
Algorithm: dynamic programming with the l2 cost.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(2)

write_data = False
outfile = "xdynp_data.txt"

n = 600
n_bkps = 3
noise_std = 1.0

signal, true_bkps = rpt.pw_constant(
    n_samples=n,
    n_bkps=n_bkps,
    noise_std=noise_std
)

if write_data:
    np.savetxt(outfile, signal, fmt="%.6f")

min_size = max(5, n // 20)
algo = rpt.Dynp(model="l2", min_size=min_size).fit(signal)
bkps_est = algo.predict(n_bkps=n_bkps)

print("n =", n)
print("true bkps =", true_bkps)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

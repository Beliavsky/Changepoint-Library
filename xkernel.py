"""
Detect changepoints in a nonlinear or non-piecewise-constant signal.
Algorithm: kernel changepoint detection with an RBF kernel.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(6)

write_data = False
outfile = "xkernel_data.txt"

n = 15000
n_bkps = 3
noise_std = 0.2

signal, true_bkps = rpt.pw_wavy(
    n_samples=n,
    n_bkps=n_bkps,
    noise_std=noise_std
)

if write_data:
    np.savetxt(outfile, signal, fmt="%.6f")

algo = rpt.KernelCPD(kernel="rbf").fit(signal)
bkps_est = algo.predict(n_bkps=n_bkps)

print("n =", n)
print("true bkps =", true_bkps)
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

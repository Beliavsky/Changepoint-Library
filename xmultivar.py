"""
Detect common changepoints in the mean structure of a multivariate series.
Algorithm: PELT with the l2 cost applied to a multivariate signal.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(5)

write_data = False
outfile = "xmultivar_data.txt"

n = 9000

frac_bkps = np.array([1/3, 2/3])
bkps_true = np.round(n * frac_bkps).astype(int)
bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))

lengths = np.diff(np.r_[0, bkps_true, n])

means_x1 = [0.0, 2.0, -1.0]
means_x2 = [1.0, 1.0, 4.0]
sd_x1 = [1.0, 1.0, 1.0]
sd_x2 = [1.0, 1.0, 1.0]

x1 = np.concatenate([
    np.random.normal(means_x1[i], sd_x1[i], lengths[i])
    for i in range(len(lengths))
])

x2 = np.concatenate([
    np.random.normal(means_x2[i], sd_x2[i], lengths[i])
    for i in range(len(lengths))
])

signal = np.column_stack([x1, x2])

if write_data:
    np.savetxt(outfile, signal, fmt="%.6f")

min_size = max(5, n // 20)
algo = rpt.Pelt(model="l2", min_size=min_size).fit(signal)
bkps_est = algo.predict(pen=12)

print("n =", n)
print("true bkps =", bkps_true.tolist() + [n])
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

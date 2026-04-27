"""
Detect changes in the mean level of a univariate series.
Algorithm: PELT with the l2 cost.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(1)

write_data = False
outfile = "xmean_shift_data.txt"

n = 300

frac_bkps = np.array([1/3, 7/12])
bkps_true = np.round(n * frac_bkps).astype(int)
bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))

lengths = np.diff(np.r_[0, bkps_true, n])

means = [0.0, 3.0, -1.0]
sds = [1.0, 1.0, 1.0]

y = np.concatenate([
    np.random.normal(means[i], sds[i], lengths[i])
    for i in range(len(lengths))
])

if write_data:
    np.savetxt(outfile, y, fmt="%.6f")

min_size = max(5, n // 20)
algo = rpt.Pelt(model="l2", min_size=min_size).fit(y)
bkps_est = algo.predict(pen=10)

print("n =", n)
print("true bkps =", bkps_true.tolist() + [n])
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

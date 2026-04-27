"""
Detect changepoints in the variance regime of a return series.
Algorithm: PELT with the radial-basis-function cost applied to squared returns.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(8)

write_data = False
outfile = "xvariance_regime_data.txt"

n = 5500

# breakpoint locations as fractions of n
frac_bkps = np.array([0.36, 0.64])
bkps_true = np.round(n * frac_bkps).astype(int)
bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))

lengths = np.diff(np.r_[0, bkps_true, n])

# standard deviations for the three variance regimes
sigmas = [0.005, 0.020, 0.010]

r = np.concatenate([
    np.random.normal(0.0, sigmas[i], lengths[i])
    for i in range(len(lengths))
])

x = (r ** 2).reshape(-1, 1)

if write_data:
    np.savetxt(outfile, x, fmt="%.6f")

min_size = max(5, n // 20)
algo = rpt.Pelt(model="rbf", min_size=min_size).fit(x)
bkps_est = algo.predict(pen=5)

print("n =", n)
print("true bkps =", bkps_true.tolist() + [n])
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

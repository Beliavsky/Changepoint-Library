"""
Detect topic or distribution shifts in a sequence of feature vectors,
such as sentence embeddings.
Algorithm: PELT with the radial-basis-function cost.
"""

import time
import numpy as np
import ruptures as rpt

t0 = time.perf_counter()
np.random.seed(7)

write_data = False
outfile = "xtext_like_data.txt"

n = 125
p = 8

frac_bkps = np.array([0.32, 0.60])
bkps_true = np.round(n * frac_bkps).astype(int)
bkps_true = np.unique(np.clip(bkps_true, 1, n - 1))

lengths = np.diff(np.r_[0, bkps_true, n])

means = [0.0, 2.0, -1.0]

emb = np.vstack([
    np.random.normal(means[i], 1.0, size=(lengths[i], p))
    for i in range(len(lengths))
])

if write_data:
    np.savetxt(outfile, emb, fmt="%.6f")

min_size = max(5, n // 20)
algo = rpt.Pelt(model="rbf", min_size=min_size).fit(emb)
bkps_est = algo.predict(pen=15)

print("n =", n)
print("p =", p)
print("true bkps =", bkps_true.tolist() + [n])
print("estimated bkps =", bkps_est)
print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")

"""
Evaluate changepoint estimates with ruptures metrics on simulated data.
"""

from __future__ import annotations

import time

import numpy as np
import ruptures as rpt
from ruptures.metrics import hausdorff, precision_recall, randindex

t0 = time.perf_counter()
np.random.seed(12)

n = 1200
n_bkps = 3
noise_std = 1.5
min_size = 40
margin = 40

signal, true_bkps = rpt.pw_constant(
    n_samples=n,
    n_bkps=n_bkps,
    noise_std=noise_std,
)

methods = [
    ("binseg", rpt.Binseg(model="l2", min_size=min_size)),
    ("bottomup", rpt.BottomUp(model="l2", min_size=min_size)),
    ("window", rpt.Window(width=max(80, n // 10), model="l2")),
    ("dynp", rpt.Dynp(model="l2", min_size=min_size, jump=1)),
]

print("n =", n)
print("true bkps =", true_bkps)
print("n_bkps =", n_bkps)
print("min_size =", min_size)
print("margin =", margin)
print()

for name, algo in methods:
    t1 = time.perf_counter()
    bkps = algo.fit(signal).predict(n_bkps=n_bkps)
    elapsed = time.perf_counter() - t1
    precision, recall = precision_recall(true_bkps, bkps, margin=margin)
    hd = hausdorff(true_bkps, bkps)
    ri = randindex(true_bkps, bkps)

    print("method =", name)
    print("estimated bkps =", bkps)
    print(f"precision = {precision:.3f}")
    print(f"recall = {recall:.3f}")
    print(f"hausdorff = {hd:.3f}")
    print(f"randindex = {ri:.6f}")
    print(f"elapsed seconds = {elapsed:.3f}")
    print()

print(f"wall time elapsed (s) = {time.perf_counter() - t0:.3f}")

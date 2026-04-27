#!/usr/bin/env python
"""Helpers for the local `bocd` package workflow."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


BOCD_REPO = Path(r"c:\python\public_domain\github\bocd")
if str(BOCD_REPO) not in sys.path:
    sys.path.insert(0, str(BOCD_REPO))

import bocd  # type: ignore  # noqa: E402


def simulate_piecewise_normal(
    n: int,
    regime_starts: list[int],
    means: list[float],
    sds: list[float],
    seed: int,
) -> tuple[np.ndarray, list[int]]:
    """Simulate a univariate piecewise-Gaussian series."""
    rng = np.random.default_rng(seed)
    x = np.empty(n, dtype=float)
    true_cps: list[int] = []
    for i, start in enumerate(regime_starts):
        i1 = start
        i2 = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        x[i1:i2] = means[i] + sds[i] * rng.normal(size=i2 - i1)
        if i + 1 < len(regime_starts):
            true_cps.append(i2)
    return x, true_cps


def write_series_file(path: str | Path, x: np.ndarray, true_cps: list[int]) -> None:
    """Write a one-column series file with changepoint metadata."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# bocd StudentT comparison data\n")
        fh.write("# true_bkps = " + " ".join(str(cp) for cp in true_cps) + "\n")
        for value in x:
            fh.write(f"{value:.12f}\n")


def read_series_file(path: str | Path) -> np.ndarray:
    """Read a one-column numeric series, skipping blank and comment lines."""
    values: list[float] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            values.append(float(line))
    return np.asarray(values, dtype=float)


def map_changepoints_bocd(runlen_argmax: np.ndarray) -> list[int]:
    """Reconstruct changepoints from the per-time MAP run length sequence."""
    cps: list[int] = []
    t = runlen_argmax.size
    while t > 0:
        r = int(runlen_argmax[t - 1])
        if r <= 0:
            t -= 1
        else:
            cp = t - r
            if cp > 0:
                cps.append(cp)
            t = cp
    cps.reverse()
    deduped: list[int] = []
    for cp in cps:
        if not deduped or deduped[-1] != cp:
            deduped.append(cp)
    return deduped


def run_bocd_student_t(
    x: np.ndarray,
    lam: float,
    mu: float = 0.0,
    kappa: float = 1.0,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> dict[str, object]:
    """Run the local bocd package with ConstantHazard and StudentT."""
    bc = bocd.BayesianOnlineChangePointDetection(
        bocd.ConstantHazard(lam),
        bocd.StudentT(mu=mu, kappa=kappa, alpha=alpha, beta=beta),
    )
    runlen_argmax = np.zeros(x.size, dtype=int)
    for i, xi in enumerate(x):
        bc.update(float(xi))
        runlen_argmax[i] = int(np.asarray(bc.rt).ravel()[0])
    map_cps = map_changepoints_bocd(runlen_argmax)
    return {
        "map_cps": map_cps,
        "runlen_argmax": runlen_argmax,
        "final_rt": int(runlen_argmax[-1]) if x.size else 0,
    }

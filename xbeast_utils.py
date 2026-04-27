#!/usr/bin/env python
"""Utilities for a small trend-only BEAST-style workflow."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def simulate_piecewise_linear(
    n: int,
    regime_starts: list[int],
    intercepts: list[float],
    slopes: list[float],
    sigma: float,
    seed: int,
) -> tuple[np.ndarray, list[int]]:
    """Simulate a univariate piecewise-linear series with Gaussian noise."""
    rng = np.random.default_rng(seed)
    y = np.empty(n, dtype=float)
    true_cps: list[int] = []
    for i, start in enumerate(regime_starts):
        i1 = start
        i2 = regime_starts[i + 1] if i + 1 < len(regime_starts) else n
        x = np.arange(i1 + 1, i2 + 1, dtype=float)
        y[i1:i2] = intercepts[i] + slopes[i] * x + sigma * rng.normal(size=i2 - i1)
        if i + 1 < len(regime_starts):
            true_cps.append(i2)
    return y, true_cps


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


def read_time_series_file(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column time series file, skipping blank and comment lines."""
    t_values: list[float] = []
    y_values: list[float] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            t_str, y_str = line.split()[:2]
            t_values.append(float(t_str))
            y_values.append(float(y_str))
    return np.asarray(t_values, dtype=float), np.asarray(y_values, dtype=float)


def read_matrix_file(path: str | Path) -> np.ndarray:
    """Read a whitespace-delimited numeric matrix, skipping blank and comment lines."""
    rows: list[list[float]] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(tok) for tok in line.split()])
    return np.asarray(rows, dtype=float)


def write_series_file(path: str | Path, x: np.ndarray, true_cps: list[int]) -> None:
    """Write a one-column series file with changepoint metadata."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# trend-only BEAST-style comparison data\n")
        fh.write("# true_bkps = " + " ".join(str(cp) for cp in true_cps) + "\n")
        for value in x:
            fh.write(f"{value:.12f}\n")


def write_time_series_file(path: str | Path, t: np.ndarray, x: np.ndarray, true_cps: list[int]) -> None:
    """Write a two-column irregular time series file with changepoint metadata."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# trend-only BEAST-style irregular comparison data\n")
        fh.write("# true_bkps = " + " ".join(str(cp) for cp in true_cps) + "\n")
        for ti, xi in zip(t, x):
            fh.write(f"{ti:.12f} {xi:.12f}\n")


def write_matrix_file(path: str | Path, x: np.ndarray, true_cps: list[int]) -> None:
    """Write a whitespace-delimited matrix file with changepoint metadata."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# beast123-style regular multi-series comparison data\n")
        fh.write("# true_bkps = " + " ".join(str(cp) for cp in true_cps) + "\n")
        for row in x:
            fh.write(" ".join(f"{value:.12f}" for value in row) + "\n")


def segment_linear_rss(prefix: dict[str, np.ndarray], a: int, b: int) -> float:
    """Return the OLS RSS for a 1-based inclusive segment [a, b]."""
    nseg = float(b - a + 1)
    sumx = prefix["sx"][b] - prefix["sx"][a - 1]
    sumx2 = prefix["sx2"][b] - prefix["sx2"][a - 1]
    sumy = prefix["sy"][b] - prefix["sy"][a - 1]
    sumy2 = prefix["sy2"][b] - prefix["sy2"][a - 1]
    sumxy = prefix["sxy"][b] - prefix["sxy"][a - 1]
    denom = nseg * sumx2 - sumx * sumx
    if abs(denom) <= 1.0e-12:
        beta1 = 0.0
        beta0 = sumy / nseg
    else:
        beta1 = (nseg * sumxy - sumx * sumy) / denom
        beta0 = (sumy - beta1 * sumx) / nseg
    rss = sumy2 - beta0 * sumy - beta1 * sumxy
    return max(float(rss), 0.0)


def fit_piecewise_linear(y: np.ndarray, cps: list[int]) -> np.ndarray:
    """Fit independent linear trends on each segment defined by changepoints."""
    n = y.size
    fitted = np.empty(n, dtype=float)
    bounds = [0] + cps + [n]
    for left, right in zip(bounds[:-1], bounds[1:]):
        x = np.arange(left + 1, right + 1, dtype=float)
        seg = y[left:right]
        xm = float(np.mean(x))
        ym = float(np.mean(seg))
        denom = float(np.sum((x - xm) ** 2))
        if denom <= 1.0e-12:
            beta1 = 0.0
            beta0 = ym
        else:
            beta1 = float(np.sum((x - xm) * (seg - ym)) / denom)
            beta0 = ym - beta1 * xm
        fitted[left:right] = beta0 + beta1 * x
    return fitted


def solve_beast_trend_only(y: np.ndarray, max_cp: int, min_seg_len: int) -> dict[str, np.ndarray | list[int] | int]:
    """Solve a small trend-only BEAST-style BMA over piecewise-linear models."""
    n = y.size
    x = np.arange(1, n + 1, dtype=float)
    prefix = {
        "sx": np.concatenate(([0.0], np.cumsum(x))),
        "sx2": np.concatenate(([0.0], np.cumsum(x * x))),
        "sy": np.concatenate(([0.0], np.cumsum(y))),
        "sy2": np.concatenate(([0.0], np.cumsum(y * y))),
        "sxy": np.concatenate(([0.0], np.cumsum(x * y))),
    }
    max_m = max_cp + 1
    dp_cost = np.full((n + 1, max_m + 1), np.inf, dtype=float)
    parent = np.zeros((n + 1, max_m + 1), dtype=int)
    for i in range(min_seg_len, n + 1):
        dp_cost[i, 1] = segment_linear_rss(prefix, 1, i)
    for m in range(2, max_m + 1):
        for i in range(m * min_seg_len, n + 1):
            best_cost = math.inf
            best_parent = 0
            for k in range((m - 1) * min_seg_len, i - min_seg_len + 1):
                cand = dp_cost[k, m - 1] + segment_linear_rss(prefix, k + 1, i)
                if cand < best_cost:
                    best_cost = cand
                    best_parent = k
            dp_cost[i, m] = best_cost
            parent[i, m] = best_parent

    model_cps: list[list[int]] = []
    fitted_models = np.zeros((max_m, n), dtype=float)
    bic = np.full(max_m, np.inf, dtype=float)
    for m in range(1, max_m + 1):
        if not np.isfinite(dp_cost[n, m]):
            model_cps.append([])
            continue
        cps: list[int] = []
        t = n
        for idx in range(m, 1, -1):
            cps.append(int(parent[t, idx]))
            t = int(parent[t, idx])
        cps.reverse()
        model_cps.append(cps)
        fitted_models[m - 1] = fit_piecewise_linear(y, cps)
        rss = max(float(np.sum((y - fitted_models[m - 1]) ** 2)), 1.0e-12)
        bic[m - 1] = n * math.log(rss / n) + (2 * m) * math.log(n)

    finite = np.isfinite(bic)
    shifted = bic[finite] - np.min(bic[finite])
    weights = np.zeros_like(bic)
    weights[finite] = np.exp(-0.5 * shifted)
    weights /= np.sum(weights)

    best_model = int(np.argmax(weights))
    cp_prob = np.zeros(n, dtype=float)
    fitted_mean = np.sum(weights[:, None] * fitted_models, axis=0)
    for m, cps in enumerate(model_cps):
        for cp in cps:
            cp_prob[cp - 1] += weights[m]

    return {
        "best_cps": model_cps[best_model],
        "ncp_best": len(model_cps[best_model]),
        "cp_prob": cp_prob,
        "fitted_mean": fitted_mean,
        "model_weights": weights,
        "bic": bic,
    }


def top_cp_probabilities(cp_prob: np.ndarray, n_top: int = 5) -> list[tuple[int, float]]:
    """Return the largest changepoint probabilities as 1-based indices."""
    work = np.array(cp_prob, copy=True)
    out: list[tuple[int, float]] = []
    for _ in range(n_top):
        pos = int(np.argmax(work))
        if work[pos] <= 0.0:
            break
        out.append((int(pos + 1), float(work[pos])))
        work[pos] = -1.0
    out.sort(key=lambda item: item[0])
    return out


def simulate_piecewise_linear_harmonic(
    n: int,
    regime_starts: list[int],
    intercepts: list[float],
    slopes: list[float],
    sigma: float,
    period: float,
    amp_sin: float,
    amp_cos: float,
    seed: int,
) -> tuple[np.ndarray, list[int]]:
    """Simulate a piecewise-linear series with a global harmonic seasonal term."""
    y, true_cps = simulate_piecewise_linear(n, regime_starts, intercepts, slopes, sigma, seed)
    x = np.arange(1, n + 1, dtype=float)
    omega = 2.0 * math.pi / period
    y = y + amp_sin * np.sin(omega * x) + amp_cos * np.cos(omega * x)
    return y, true_cps


def simulate_irregular_time_grid(
    n: int,
    start_time: float,
    dt_mean: float,
    dt_jitter: float,
    seed: int,
) -> np.ndarray:
    """Simulate a strictly increasing irregular time grid."""
    rng = np.random.default_rng(seed)
    t = np.empty(n, dtype=float)
    t[0] = start_time
    increments = dt_mean + dt_jitter * (2.0 * rng.random(n - 1) - 1.0)
    increments = np.maximum(increments, 0.1 * dt_mean)
    t[1:] = start_time + np.cumsum(increments)
    return t


def simulate_piecewise_linear_irregular(
    t: np.ndarray,
    regime_starts: list[int],
    intercepts: list[float],
    slopes: list[float],
    sigma: float,
    seed: int,
) -> tuple[np.ndarray, list[int]]:
    """Simulate an irregular-time piecewise-linear series."""
    rng = np.random.default_rng(seed)
    y = np.empty(t.size, dtype=float)
    true_cps: list[int] = []
    for i, start in enumerate(regime_starts):
        i1 = start
        i2 = regime_starts[i + 1] if i + 1 < len(regime_starts) else t.size
        y[i1:i2] = intercepts[i] + slopes[i] * t[i1:i2] + sigma * rng.normal(size=i2 - i1)
        if i + 1 < len(regime_starts):
            true_cps.append(i2)
    return y, true_cps


def fit_piecewise_linear_irregular(t: np.ndarray, y: np.ndarray, cps: list[int]) -> np.ndarray:
    """Fit independent linear trends on each irregular-time segment."""
    fitted = np.empty(y.size, dtype=float)
    bounds = [0] + cps + [y.size]
    for left, right in zip(bounds[:-1], bounds[1:]):
        xseg = t[left:right]
        yseg = y[left:right]
        xm = float(np.mean(xseg))
        ym = float(np.mean(yseg))
        denom = float(np.sum((xseg - xm) ** 2))
        if denom <= 1.0e-12:
            beta1 = 0.0
            beta0 = ym
        else:
            beta1 = float(np.sum((xseg - xm) * (yseg - ym)) / denom)
            beta0 = ym - beta1 * xm
        fitted[left:right] = beta0 + beta1 * xseg
    return fitted


def segment_linear_rss_irregular(prefix: dict[str, np.ndarray], a: int, b: int) -> float:
    """Return the OLS RSS for an irregular-time 1-based inclusive segment [a, b]."""
    nseg = float(b - a + 1)
    sumx = prefix["sx"][b] - prefix["sx"][a - 1]
    sumx2 = prefix["sx2"][b] - prefix["sx2"][a - 1]
    sumy = prefix["sy"][b] - prefix["sy"][a - 1]
    sumy2 = prefix["sy2"][b] - prefix["sy2"][a - 1]
    sumxy = prefix["sxy"][b] - prefix["sxy"][a - 1]
    denom = nseg * sumx2 - sumx * sumx
    if abs(denom) <= 1.0e-12:
        beta1 = 0.0
        beta0 = sumy / nseg
    else:
        beta1 = (nseg * sumxy - sumx * sumy) / denom
        beta0 = (sumy - beta1 * sumx) / nseg
    rss = sumy2 - beta0 * sumy - beta1 * sumxy
    return max(float(rss), 0.0)


def solve_beast_irreg_trend_only(t: np.ndarray, y: np.ndarray, max_cp: int, min_seg_len: int) -> dict[str, np.ndarray | list[int] | int]:
    """Solve an irregular-time trend-only BEAST-style BMA."""
    n = y.size
    prefix = {
        "sx": np.concatenate(([0.0], np.cumsum(t))),
        "sx2": np.concatenate(([0.0], np.cumsum(t * t))),
        "sy": np.concatenate(([0.0], np.cumsum(y))),
        "sy2": np.concatenate(([0.0], np.cumsum(y * y))),
        "sxy": np.concatenate(([0.0], np.cumsum(t * y))),
    }
    max_m = max_cp + 1
    dp_cost = np.full((n + 1, max_m + 1), np.inf, dtype=float)
    parent = np.zeros((n + 1, max_m + 1), dtype=int)
    for i in range(min_seg_len, n + 1):
        dp_cost[i, 1] = segment_linear_rss_irregular(prefix, 1, i)
    for m in range(2, max_m + 1):
        for i in range(m * min_seg_len, n + 1):
            best_cost = math.inf
            best_parent = 0
            for k in range((m - 1) * min_seg_len, i - min_seg_len + 1):
                cand = dp_cost[k, m - 1] + segment_linear_rss_irregular(prefix, k + 1, i)
                if cand < best_cost:
                    best_cost = cand
                    best_parent = k
            dp_cost[i, m] = best_cost
            parent[i, m] = best_parent

    model_cps: list[list[int]] = []
    fitted_models = np.zeros((max_m, n), dtype=float)
    bic = np.full(max_m, np.inf, dtype=float)
    for m in range(1, max_m + 1):
        if not np.isfinite(dp_cost[n, m]):
            model_cps.append([])
            continue
        cps: list[int] = []
        idx = n
        for col in range(m, 1, -1):
            cps.append(int(parent[idx, col]))
            idx = int(parent[idx, col])
        cps.reverse()
        model_cps.append(cps)
        fitted_models[m - 1] = fit_piecewise_linear_irregular(t, y, cps)
        rss = max(float(np.sum((y - fitted_models[m - 1]) ** 2)), 1.0e-12)
        bic[m - 1] = n * math.log(rss / n) + (2 * m) * math.log(n)

    finite = np.isfinite(bic)
    shifted = bic[finite] - np.min(bic[finite])
    weights = np.zeros_like(bic)
    weights[finite] = np.exp(-0.5 * shifted)
    weights /= np.sum(weights)

    best_model = int(np.argmax(weights))
    cp_prob = np.zeros(n, dtype=float)
    fitted_mean = np.sum(weights[:, None] * fitted_models, axis=0)
    for m, cps in enumerate(model_cps):
        for cp in cps:
            cp_prob[cp - 1] += weights[m]

    return {
        "best_cps": model_cps[best_model],
        "ncp_best": len(model_cps[best_model]),
        "cp_prob": cp_prob,
        "fitted_mean": fitted_mean,
        "model_weights": weights,
        "bic": bic,
    }


def solve_beast123_regular(x: np.ndarray, max_cp: int, min_seg_len: int) -> dict[str, object]:
    """Apply the regular trend-only BEAST-style solver column-wise to multiple series."""
    n, p = x.shape
    best_cps: list[list[int]] = []
    model_weights = np.zeros((p, max_cp + 1), dtype=float)
    cp_prob = np.zeros((n, p), dtype=float)
    fitted_mean = np.zeros((n, p), dtype=float)
    for j in range(p):
        out = solve_beast_trend_only(x[:, j], max_cp=max_cp, min_seg_len=min_seg_len)
        best_cps.append(list(out["best_cps"]))
        model_weights[j, :] = np.asarray(out["model_weights"])
        cp_prob[:, j] = np.asarray(out["cp_prob"])
        fitted_mean[:, j] = np.asarray(out["fitted_mean"])
    mean_cp_prob = np.mean(cp_prob, axis=1)
    return {
        "best_cps": best_cps,
        "model_weights": model_weights,
        "cp_prob": cp_prob,
        "mean_cp_prob": mean_cp_prob,
        "fitted_mean": fitted_mean,
    }


def fit_harmonic_component(y: np.ndarray, period: float) -> tuple[np.ndarray, np.ndarray]:
    """Fit a fixed-period harmonic component a*sin + b*cos by least squares."""
    x = np.arange(1, y.size + 1, dtype=float)
    omega = 2.0 * math.pi / period
    design = np.column_stack((np.sin(omega * x), np.cos(omega * x)))
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    seasonal = design @ beta
    return seasonal, beta


def solve_beast_harmonic(y: np.ndarray, period: float, max_cp: int, min_seg_len: int) -> dict[str, np.ndarray | list[int] | int]:
    """Solve a harmonic-seasonal plus trend-changepoint BEAST-style approximation."""
    seasonal, beta = fit_harmonic_component(y, period)
    deseasonalized = y - seasonal
    out = solve_beast_trend_only(deseasonalized, max_cp=max_cp, min_seg_len=min_seg_len)
    out["seasonal_fit"] = seasonal
    out["beta_season"] = beta
    out["deseasonalized"] = deseasonalized
    out["fitted_mean"] = np.asarray(out["fitted_mean"]) + seasonal
    return out

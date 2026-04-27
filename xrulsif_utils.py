"""
Shared helpers for deterministic Gaussian-kernel RuLSIF comparison scripts.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from xsst_utils import read_series_file, write_series_file


def simulate_rulsif_signal(
    n: int,
    true_cp: int,
    means: tuple[float, float],
    sds: tuple[float, float],
    seed: int,
) -> np.ndarray:
    """
    Simulate a univariate Gaussian series with one changepoint.
    """
    rng = np.random.default_rng(seed)
    signal = np.empty(n, dtype=float)
    signal[:true_cp] = rng.normal(means[0], sds[0], size=true_cp)
    signal[true_cp:] = rng.normal(means[1], sds[1], size=n - true_cp)
    return signal


def compile_hankel_1d(signal: np.ndarray, end_index: int, window_length: int, n_cols: int) -> np.ndarray:
    """
    Build the same slow Hankel matrix as changepoynt.utils.linalg.compile_hankel with lag=1.
    """
    hankel = np.empty((window_length, n_cols), dtype=float)
    for cx in range(n_cols):
        hankel[:, -cx - 1] = signal[end_index - window_length - cx : end_index - cx]
    return hankel


def gaussian_kernel_matrix(samples: np.ndarray, centers: np.ndarray, sigma: float) -> np.ndarray:
    """
    Compute the Gaussian kernel matrix with center-major row layout.
    """
    sample_sq = np.sum(samples * samples, axis=0)
    center_sq = np.sum(centers * centers, axis=0)
    sqdist = center_sq[:, None] + sample_sq[None, :] - 2.0 * centers.T @ samples
    return np.exp(-sqdist / (2.0 * sigma * sigma))


def rulsif_pair_score(
    reference_samples: np.ndarray,
    test_samples: np.ndarray,
    alpha: float,
    sigma: float,
    lambda_reg: float,
) -> float:
    """
    Deterministic RuLSIF PE-divergence score for one reference/test pair.
    """
    all_samples = np.concatenate((reference_samples, test_samples), axis=1)
    std = np.std(all_samples, axis=1) + np.finfo(float).eps
    std[std < 1e-12] = 1.0
    ref_scaled = reference_samples / std[:, None]
    test_scaled = test_samples / std[:, None]

    centers = ref_scaled
    k_ref = gaussian_kernel_matrix(ref_scaled, centers, sigma)
    k_test = gaussian_kernel_matrix(test_scaled, centers, sigma)

    n_ref = ref_scaled.shape[1]
    n_test = test_scaled.shape[1]
    h_hat = np.mean(k_ref, axis=1)
    h_mat = alpha / n_ref * (k_ref @ k_ref.T) + (1.0 - alpha) / n_test * (k_test @ k_test.T)
    h_mat = h_mat + lambda_reg * np.eye(h_mat.shape[0])
    theta = np.linalg.solve(h_mat, h_hat)
    ratios_test = np.maximum(theta @ k_test, 0.0)
    return float(0.5 * np.mean(ratios_test) - 0.5)


def fit_rulsif_gaussian(
    signal: np.ndarray,
    window_length: int,
    n_windows: int,
    lag: int,
    alpha: float,
    sigma: float,
    lambda_reg: float,
    scoring_step: int = 1,
    symmetric: bool = True,
) -> tuple[np.ndarray, int, float]:
    """
    Run the deterministic Gaussian-kernel RuLSIF score path.
    """
    signal = np.asarray(signal, dtype=float)
    start_idx = window_length + n_windows + lag
    score = np.zeros_like(signal)
    offset = n_windows

    for idx in range(start_idx, signal.shape[0], scoring_step):
        hankel = compile_hankel_1d(signal, idx, window_length, 2 * n_windows)
        local_score = rulsif_pair_score(
            hankel[:, :n_windows],
            hankel[:, n_windows:],
            alpha=alpha,
            sigma=sigma,
            lambda_reg=lambda_reg,
        )
        left_idx = max(idx - offset - scoring_step // 2, 0)
        right_idx = min(idx - offset + (scoring_step + 1) // 2, signal.shape[0])
        score[left_idx:right_idx] = local_score

    if symmetric:
        rev_score, _, _ = fit_rulsif_gaussian(
            signal[::-1],
            window_length=window_length,
            n_windows=n_windows,
            lag=lag,
            alpha=alpha,
            sigma=sigma,
            lambda_reg=lambda_reg,
            scoring_step=scoring_step,
            symmetric=False,
        )
        score = score + rev_score[::-1]

    cp = int(np.argmax(score))
    return score, cp, float(score[cp])

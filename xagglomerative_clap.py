"""
Simulate a repeated-state univariate series and run AgglomerativeCLaP-style merging.

This mirrors claspy.state_detection.AgglomerativeCLaPDetection, but runs on top of
the deterministic CLaP-style classifier used in this repo so Python and Fortran can
be compared exactly.
"""

from __future__ import annotations

import time

import numpy as np

from xclap import fit_clap_centroid, macro_f1_labels
from xclasp import simulate_clasp_signal


def create_state_labels(true_bkps: list[int], segment_labels: list[int], n: int) -> np.ndarray:
    """
    Expand per-segment labels into pointwise state labels.
    """
    state_labels = np.empty(n, dtype=int)
    start = 0
    for i, stop in enumerate(true_bkps + [n]):
        state_labels[start:stop] = segment_labels[i]
        start = stop
    return state_labels


def random_f1_labels(y_true: np.ndarray) -> float:
    """
    Expected macro-F1 score under random predictions with matching class frequencies.
    """
    labels = np.unique(y_true)
    score = 0.0
    n = y_true.shape[0]
    for label in labels:
        pos = np.sum(y_true == label)
        neg = n - pos
        tp = pos * pos / n
        fn = pos * neg / n
        fp = neg * pos / n
        pre = tp / (tp + fp)
        rec = tp / (tp + fn)
        if pre + rec > 0:
            score += 2.0 * pre * rec / (pre + rec)
    return score / labels.shape[0]


def classification_gain(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Gain over the random-label macro-F1 baseline.
    """
    return macro_f1_labels(y_true, y_pred) - random_f1_labels(y_true)


def collapse_segment_process(
    true_bkps: list[int], segment_labels: list[int]
) -> tuple[list[int], list[int]]:
    """
    Collapse consecutive equal segment labels into sparse state transitions.
    """
    labels_out = [segment_labels[0]]
    cps_out: list[int] = []
    for i in range(1, len(segment_labels)):
        if segment_labels[i] != segment_labels[i - 1]:
            cps_out.append(true_bkps[i - 1])
            labels_out.append(segment_labels[i])
    return cps_out, labels_out


def fit_agglomerative_clap(
    signal: np.ndarray,
    true_bkps: list[int],
    window_size: int,
    n_splits: int,
    sample_cap: int = 1000,
) -> tuple[list[int], float, np.ndarray, np.ndarray]:
    """
    Run deterministic agglomerative merging of initially unique segment labels.
    """
    labels = np.arange(1, len(true_bkps) + 2, dtype=int)
    ignore_pairs: set[tuple[int, int]] = set()

    while True:
        state_labels = create_state_labels(true_bkps, labels.tolist(), signal.shape[0])
        y_true, y_pred, _ = fit_clap_centroid(
            signal, state_labels, window_size, n_splits, sample_cap
        )
        current_gain = classification_gain(y_true, y_pred)
        unique_labels = np.unique(labels)
        if unique_labels.shape[0] <= 1:
            break

        conf = np.zeros((unique_labels.shape[0], unique_labels.shape[0]), dtype=int)
        idx = {label: i for i, label in enumerate(unique_labels.tolist())}
        for yt, yp in zip(y_true.tolist(), y_pred.tolist(), strict=True):
            conf[idx[yt], idx[yp]] += 1

        conf_loss = np.zeros(unique_labels.shape[0], dtype=float)
        conf_index = np.zeros(unique_labels.shape[0], dtype=int)
        for i, row in enumerate(conf):
            tmp = row.copy()
            tmp[i] = 0
            conf_index[i] = int(np.argmax(tmp))
            conf_loss[i] = float(np.max(tmp) / np.sum(row)) if np.sum(row) > 0 else 0.0

        merged = False
        for order_idx in np.argsort(conf_loss)[::-1]:
            label1 = int(unique_labels[order_idx])
            label2 = int(unique_labels[conf_index[order_idx]])
            if label1 == label2:
                continue

            pair = (min(label1, label2), max(label1, label2))
            if pair in ignore_pairs:
                continue

            y_true_m = y_true.copy()
            y_pred_m = y_pred.copy()
            y_true_m[y_true_m == label2] = label1
            y_pred_m[y_pred_m == label2] = label1
            new_gain = classification_gain(y_true_m, y_pred_m)
            if current_gain > new_gain:
                ignore_pairs.add(pair)
                continue

            if label2 > label1:
                label1, label2 = label2, label1
            labels[labels == label2] = label1
            merged = True
            break

        if not merged:
            break

    mapped = {label: i + 1 for i, label in enumerate(np.unique(labels).tolist())}
    final_segment_labels = [mapped[int(label)] for label in labels.tolist()]
    state_labels = create_state_labels(true_bkps, final_segment_labels, signal.shape[0])
    y_true, y_pred, _ = fit_clap_centroid(
        signal, state_labels, window_size, n_splits, sample_cap
    )
    gain = classification_gain(y_true, y_pred)
    return final_segment_labels, gain, y_true, y_pred


def main() -> None:
    t0 = time.perf_counter()
    n = 480
    true_bkps = [80, 160, 240, 320, 400]
    true_states = [1, 2, 3, 1, 2, 3]
    means = [0.0, 2.0, -2.0, 0.0, 2.0, -2.0]
    sds = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    signal = simulate_clasp_signal(n, [0] + true_bkps, means, sds, 101)

    window_size = 20
    n_splits = 5
    sample_cap = 1000
    segment_labels, gain, y_true, y_pred = fit_agglomerative_clap(
        signal, true_bkps, window_size, n_splits, sample_cap
    )
    cps_out, labels_out = collapse_segment_process(true_bkps, segment_labels)

    print("n =", n)
    print("true_bkps =", true_bkps)
    print("true states =", true_states)
    print("window_size =", window_size)
    print("n_splits =", n_splits)
    print("sample_cap =", sample_cap)
    print("segment labels =", segment_labels)
    print("collapsed bkps =", cps_out)
    print("collapsed labels =", labels_out)
    print("y_true =", y_true.tolist())
    print("y_pred =", y_pred.tolist())
    print(f"classification gain = {gain:.6f}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()

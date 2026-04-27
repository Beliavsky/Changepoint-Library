#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xbinseg_data.txt"
n_bkps <- 4L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
min_size <- max(5L, n %/% 25L)

prefix <- build_prefix_1d(signal)
cs <- prefix$cs
css <- prefix$css

segments <- matrix(c(1L, n), ncol = 2L)
bkps <- integer(0)

for (iter in seq_len(n_bkps)) {
  best_gain <- -Inf
  best_seg <- NA_integer_
  best_k <- NA_integer_

  for (iseg in seq_len(nrow(segments))) {
    s <- segments[iseg, 1L]
    e <- segments[iseg, 2L]
    parent_cost <- segment_sse_1d(cs, css, s, e, min_size)
    if (!is.finite(parent_cost)) next

    for (k in s:(e - 1L)) {
      left_cost <- segment_sse_1d(cs, css, s, k, min_size)
      if (!is.finite(left_cost)) next
      right_cost <- segment_sse_1d(cs, css, k + 1L, e, min_size)
      if (!is.finite(right_cost)) next

      gain <- parent_cost - left_cost - right_cost
      if (gain > best_gain) {
        best_gain <- gain
        best_seg <- iseg
        best_k <- k
      }
    }
  }

  if (is.na(best_seg)) break
  bkps <- c(bkps, best_k)

  old_s <- segments[best_seg, 1L]
  old_e <- segments[best_seg, 2L]
  left_seg <- matrix(c(old_s, best_k), ncol = 2L)
  right_seg <- matrix(c(best_k + 1L, old_e), ncol = 2L)

  if (nrow(segments) == 1L) {
    segments <- rbind(left_seg, right_seg)
  } else if (best_seg == 1L) {
    segments <- rbind(left_seg, right_seg, segments[-1L, , drop = FALSE])
  } else if (best_seg == nrow(segments)) {
    segments <- rbind(segments[-best_seg, , drop = FALSE], left_seg, right_seg)
  } else {
    segments <- rbind(
      segments[seq_len(best_seg - 1L), , drop = FALSE],
      left_seg,
      right_seg,
      segments[(best_seg + 1L):nrow(segments), , drop = FALSE]
    )
  }
}

bkps <- sort(bkps)
bkps_out <- c(bkps, n)
elapsed <- proc.time()[["elapsed"]] - t0

print_result_1d(
  data_file,
  n,
  list(
    list(name = "n_bkps", value = n_bkps),
    list(name = "min_size", value = min_size)
  ),
  bkps_out,
  elapsed
)

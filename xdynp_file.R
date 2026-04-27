#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdynp_data.txt"
n_bkps <- 3L
max_m <- n_bkps + 1L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
min_size <- max(5L, n %/% 20L)

prefix <- build_prefix_1d(signal)
cs <- prefix$cs
css <- prefix$css

dp <- matrix(Inf, nrow = n, ncol = max_m)
parent <- matrix(0L, nrow = n, ncol = max_m)

for (i in seq_len(n)) {
  dp[i, 1L] <- segment_sse_1d(cs, css, 1L, i, min_size)
}

if (max_m >= 2L) {
  for (m in 2:max_m) {
    for (i in m:n) {
      best_val <- Inf
      best_k <- 0L
      for (k in (m - 1L):(i - 1L)) {
        if (!is.finite(dp[k, m - 1L])) next
        seg_cost <- segment_sse_1d(cs, css, k + 1L, i, min_size)
        cand <- dp[k, m - 1L] + seg_cost
        if (cand < best_val) {
          best_val <- cand
          best_k <- k
        }
      }
      dp[i, m] <- best_val
      parent[i, m] <- best_k
    }
  }
}

seg_ends <- integer(max_m)
idx <- n
seg_ends[max_m] <- n
if (max_m >= 2L) {
  for (m in seq(max_m, 2L, by = -1L)) {
    seg_ends[m - 1L] <- parent[idx, m]
    idx <- seg_ends[m - 1L]
  }
}

elapsed <- proc.time()[["elapsed"]] - t0

print_result_1d(
  data_file,
  n,
  list(
    list(name = "n_bkps", value = n_bkps),
    list(name = "min_size", value = min_size)
  ),
  seg_ends,
  elapsed
)

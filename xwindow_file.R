#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xwindow_data.txt"
n_bkps <- 3L
jump <- 5L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
width <- max(20L, n %/% 10L)
width <- 2L * (width %/% 2L)
min_size <- 2L

prefix <- build_prefix_1d(signal)
cs <- prefix$cs
css <- prefix$css

inds <- seq.int(0L, n - 1L, by = jump)
keep <- inds >= width %/% 2L & inds < n - width %/% 2L
inds <- inds[keep]

if (length(inds) == 0L) {
  bkps_out <- n
} else {
  score <- numeric(length(inds))
  for (i in seq_along(inds)) {
    k <- inds[i]
    start0 <- k - width %/% 2L
    end0 <- k + width %/% 2L
    full_cost <- segment_sse_1d(cs, css, start0 + 1L, end0, 1L)
    left_cost <- segment_sse_1d(cs, css, start0 + 1L, k, 1L)
    right_cost <- segment_sse_1d(cs, css, k + 1L, end0, 1L)
    score[i] <- full_cost - left_cost - right_cost
  }

  order <- max(max(width, 2L * min_size) %/% (2L * jump), 1L)
  peaks <- find_window_peaks(score, inds, order)

  if (length(peaks$bkps) == 0L) {
    bkps_out <- n
  } else {
    ord <- order(peaks$gain)
    peak_bkps <- peaks$bkps[ord]
    chosen <- tail(peak_bkps, min(n_bkps, length(peak_bkps)))
    bkps_out <- c(sort(chosen), n)
  }
}

elapsed <- proc.time()[["elapsed"]] - t0

print_result_1d(
  data_file,
  n,
  list(
    list(name = "width", value = width),
    list(name = "n_bkps", value = n_bkps)
  ),
  bkps_out,
  elapsed
)

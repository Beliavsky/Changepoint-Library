#!/usr/bin/env Rscript

source("utils.r")
library(kerSeg)

#' Find the best kerSeg split for one inclusive segment.
segment_best_split <- function(kernel, start_idx, end_idx, min_size) {
  len <- end_idx - start_idx + 1L
  if (len < 2L * min_size) {
    return(NULL)
  }

  sub_kernel <- kernel[start_idx:end_idx, start_idx:end_idx, drop = FALSE]
  fit <- kerseg1(
    n = len,
    K = sub_kernel,
    n0 = min_size,
    n1 = len - min_size,
    pval.appr = FALSE,
    skew.corr = FALSE,
    pval.perm = FALSE
  )

  tauhat <- as.integer(fit$tauhat)
  score <- suppressWarnings(as.numeric(fit$stat$GPK$Zmax))
  if (
    !is.finite(tauhat) ||
    tauhat < min_size ||
    tauhat > len - min_size ||
    length(score) != 1L ||
    !is.finite(score)
  ) {
    return(NULL)
  }

  list(
    bkp = start_idx + tauhat - 1L,
    score = score
  )
}

#' Run recursive binary segmentation using kerSeg single-change fits.
kernel_binseg <- function(kernel, n_bkps, min_size) {
  n <- nrow(kernel)
  segments <- list(c(1L, n))
  bkps <- integer(0)

  for (step in seq_len(n_bkps)) {
    best <- NULL
    best_seg <- 0L

    for (i in seq_along(segments)) {
      seg <- segments[[i]]
      cand <- segment_best_split(kernel, seg[1L], seg[2L], min_size)
      if (is.null(cand)) next
      if (is.null(best) || cand$score > best$score) {
        best <- cand
        best_seg <- i
      }
    }

    if (is.null(best)) break

    seg <- segments[[best_seg]]
    left_seg <- c(seg[1L], best$bkp)
    right_seg <- c(best$bkp + 1L, seg[2L])
    segments <- append(segments[-best_seg], list(left_seg, right_seg))
    bkps <- sort(c(bkps, best$bkp))
  }

  c(bkps, n)
}

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xkernel_data.txt"
n_bkps <- 3L
min_size <- 2L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
x <- matrix(signal, ncol = 1L)
n <- nrow(x)

kernel <- gaussiankernel(x)
bkps <- kernel_binseg(kernel, n_bkps, min_size)
elapsed <- proc.time()[["elapsed"]] - t0

print_result_1d(
  data_file,
  n,
  list(
    list(name = "n_bkps", value = n_bkps),
    list(name = "min_size", value = min_size)
  ),
  bkps,
  elapsed
)

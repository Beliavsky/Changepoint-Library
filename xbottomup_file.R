#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xbottomup_data.txt"
n_bkps <- 4L
jump <- 1L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
min_size <- max(5L, n %/% 25L)

bkps_init <- bottomup_l2_1d(signal, n_bkps, min_size, jump)
bkps_est <- refine_bkps_l2_1d(signal, bkps_init, min_size)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("n_bkps =", n_bkps, "\n")
cat("min_size =", min_size, "\n")
cat("initial bkps =", paste(bkps_init, collapse = " "), "\n")
cat("estimated bkps =", paste(bkps_est, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

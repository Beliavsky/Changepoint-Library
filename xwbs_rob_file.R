#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xwbs_rob_data.txt"
interval_file <- sub("\\.txt$", "_intervals.txt", data_file)
delta <- 5L
K <- 1.345
tau <- 8

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
intervals <- read.table(interval_file, comment.char = "#")
raw <- changepoints::WBS.uni.rob(y, 1L, length(y), intervals[[1]], intervals[[2]], K = K, delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0

cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))

cat("file =", data_file, "\n")
cat("interval_file =", interval_file, "\n")
cat("n =", length(y), "\n")
cat("M =", nrow(intervals), "\n")
cat("K =", K, "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

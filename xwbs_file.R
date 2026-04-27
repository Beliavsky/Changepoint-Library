#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xwbs_data.txt"
interval_file <- sub("\\.txt$", "_intervals.txt", data_file)
delta <- 5L
tau <- 4

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
intervals <- read.table(interval_file, comment.char = "#")
raw <- changepoints::WBS.univar(y, 1L, length(y), intervals[[1]], intervals[[2]], delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0

if (is.null(trimmed$cpt_hat)) {
    cps_hat <- integer()
} else {
    cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))
}

cat("file =", data_file, "\n")
cat("interval_file =", interval_file, "\n")
cat("n =", length(y), "\n")
cat("M =", nrow(intervals), "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

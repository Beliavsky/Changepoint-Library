#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xwbsip_cov_data.txt"
prime_file <- sub("\\.txt$", "_prime.txt", data_file)
interval_file <- sub("\\.txt$", "_intervals.txt", data_file)
delta <- 5L

X <- t(as.matrix(read.table(data_file, comment.char = "#")))
X_prime <- t(as.matrix(read.table(prime_file, comment.char = "#")))
intervals <- read.table(interval_file, comment.char = "#")
tau <- 8

t0 <- proc.time()[["elapsed"]]
raw <- changepoints::WBSIP.cov(X, X_prime, 1L, ncol(X), intervals[[1]], intervals[[2]], delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0

if (is.null(trimmed$cpt_hat)) {
    cps_hat <- integer()
} else {
    cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))
}

cat("file =", data_file, "\n")
cat("prime_file =", prime_file, "\n")
cat("interval_file =", interval_file, "\n")
cat("n =", ncol(X), "\n")
cat("p =", nrow(X), "\n")
cat("M =", nrow(intervals), "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xbs_cov_data.txt"
tau <- 7

t0 <- proc.time()[["elapsed"]]
X <- t(as.matrix(read.table(data_file, comment.char = "#")))
raw <- changepoints::BS.cov(X, 1L, ncol(X))
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0

if (is.null(trimmed$cpt_hat)) {
    cps_hat <- integer()
} else {
    cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))
}

cat("file =", data_file, "\n")
cat("n =", ncol(X), "\n")
cat("p =", nrow(X), "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

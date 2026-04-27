#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdpdu_regression_data.txt"
lambda <- 0.5
zeta <- 20L

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
X <- dat[, -1, drop = FALSE]

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::DPDU.regression(y = y, X = X, lambda = lambda, zeta = zeta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("p =", ncol(X), "\n")
cat("lambda =", lambda, "\n")
cat("zeta =", zeta, "\n")
cat("estimated changepoints =", paste(as.numeric(ans$cpt), collapse = " "), "\n")
cat(sprintf("beta checksum = %.12f\n", sum(ans$beta_mat)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

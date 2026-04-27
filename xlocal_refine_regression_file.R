#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdp_regression_data.txt"
gamma <- 2
lambda <- 0.5
delta <- 5L
zeta <- 0.5

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
X <- dat[, -1, drop = FALSE]

t0 <- proc.time()[["elapsed"]]
init <- changepoints::DP.regression(y = y, X = X, gamma = gamma, lambda = lambda, delta = delta)$cpt
refined <- changepoints::local.refine.regression(init, y = y, X = X, zeta = zeta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("p =", ncol(X), "\n")
cat("gamma =", gamma, "\n")
cat("lambda =", lambda, "\n")
cat("delta =", delta, "\n")
cat("zeta =", zeta, "\n")
cat("initial changepoints =", paste(init, collapse = " "), "\n")
cat("refined changepoints =", paste(refined, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

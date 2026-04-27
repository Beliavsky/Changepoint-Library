#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdp_var1_data.txt"
gamma <- 0.2
lambda <- 0
delta <- 5L

DATA <- t(as.matrix(read.table(data_file, comment.char = "#")))
X_curr <- DATA[, 1:(ncol(DATA) - 1), drop = FALSE]
X_futu <- DATA[, 2:ncol(DATA), drop = FALSE]

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::DP.VAR1(X_futu = X_futu, X_curr = X_curr, gamma = gamma, lambda = lambda, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("p =", nrow(DATA), "\n")
cat("n transitions =", ncol(X_curr), "\n")
cat("gamma =", gamma, "\n")
cat("lambda =", lambda, "\n")
cat("delta =", delta, "\n")
cat("estimated changepoints =", paste(ans$cpt, collapse = " "), "\n")
cat(sprintf("partition checksum = %.12f\n", sum(ans$partition)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

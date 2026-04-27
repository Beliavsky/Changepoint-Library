#!/usr/bin/env Rscript

set.seed(0)
n <- 300L
gamma <- 5
delta <- 5L
true_cps <- c(20L, 50L, 170L)
y <- c(rep(0, 20), rep(2, 30), rep(0, 120), rep(-2, 130)) + rnorm(n, 0, 1)

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::DP.univar(y, gamma = gamma, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("gamma =", gamma, "\n")
cat("delta =", delta, "\n")
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
cat("estimated changepoints =", paste(ans$cpt, collapse = " "), "\n")
cat("n changepoints =", length(ans$cpt), "\n")
cat("partition last =", tail(ans$partition, 1), "\n")
cat(sprintf("yhat checksum = %.12f\n", sum(ans$yhat)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

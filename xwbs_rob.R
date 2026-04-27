#!/usr/bin/env Rscript

set.seed(0)
n <- 300L
delta <- 5L
K <- 1.345
tau <- 8
M <- 300L
true_cps <- c(20L, 50L, 170L)
y <- c(rep(0, 20), rep(2, 30), rep(0, 120), rep(-2, 130)) + rnorm(n, 0, 1)
outlier_idx <- c(15L, 35L, 60L, 90L, 140L, 180L, 220L, 260L, 290L)
outlier_shift <- c(15, -12, 14, -16, 13, 12, -14, 15, -13)
y[outlier_idx] <- y[outlier_idx] + outlier_shift
set.seed(1)
intervals <- changepoints::WBS.intervals(M = M, lower = 1L, upper = n)

t0 <- proc.time()[["elapsed"]]
raw <- changepoints::WBS.uni.rob(y, 1L, n, intervals$Alpha, intervals$Beta, K = K, delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0

cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))

cat("n =", n, "\n")
cat("M =", M, "\n")
cat("K =", K, "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
cat("outlier_idx =", paste(outlier_idx, collapse = " "), "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

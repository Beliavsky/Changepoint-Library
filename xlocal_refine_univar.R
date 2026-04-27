#!/usr/bin/env Rscript

set.seed(0)
n <- 300L
delta <- 5L
tau <- 4
M <- 300L
true_cps <- c(20L, 50L, 170L)
y <- c(rep(0, 20), rep(2, 30), rep(0, 120), rep(-2, 130)) + rnorm(n, 0, 1)
set.seed(1)
intervals <- changepoints::WBS.intervals(M = M, lower = 1L, upper = n)

t0 <- proc.time()[["elapsed"]]
raw <- changepoints::WBS.univar(y, 1L, n, intervals$Alpha, intervals$Beta, delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
if (is.null(trimmed$cpt_hat)) {
    cps_init <- integer()
    cps_refined <- integer()
} else {
    cps_init <- sort(as.integer(trimmed$cpt_hat[, 1]))
    cps_refined <- changepoints::local.refine.univar(cps_init, y)
}
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("M =", M, "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
cat("initial changepoints =", paste(cps_init, collapse = " "), "\n")
cat("refined changepoints =", paste(cps_refined, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

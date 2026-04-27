#!/usr/bin/env Rscript

set.seed(1234)
n <- 180L
p <- 4L
M <- 200L
delta <- 5L
A1 <- matrix(c(
    1.0, 0.65, 0.0, 0.0,
    0.65, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.15,
    0.0, 0.0, 0.15, 1.0
), nrow = p, byrow = TRUE)
A2 <- matrix(c(
    1.7, -0.15, 0.0, 0.0,
    -0.15, 0.8, 0.0, 0.0,
    0.0, 0.0, 0.7, -0.55,
    0.0, 0.0, -0.55, 1.5
), nrow = p, byrow = TRUE)
X <- cbind(
    t(MASS::mvrnorm(80, mu = rep(0, p), Sigma = A1)),
    t(MASS::mvrnorm(100, mu = rep(0, p), Sigma = A2))
)
set.seed(4321)
X_prime <- cbind(
    t(MASS::mvrnorm(80, mu = rep(0, p), Sigma = A1)),
    t(MASS::mvrnorm(100, mu = rep(0, p), Sigma = A2))
)
set.seed(202)
intervals <- changepoints::WBS.intervals(M = M, lower = 1, upper = ncol(X))
tau <- 8

t0 <- proc.time()[["elapsed"]]
raw <- changepoints::WBSIP.cov(X, X_prime, 1L, ncol(X), intervals$Alpha, intervals$Beta, delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0
if (is.null(trimmed$cpt_hat)) {
    cps_hat <- integer()
} else {
    cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))
}

cat("n =", n, "\n")
cat("p =", p, "\n")
cat("M =", M, "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

set.seed(1234)
n <- 180L
p <- 4L
tau <- 7
A1 <- matrix(c(
    1.0, 0.7, 0.0, 0.0,
    0.7, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.2,
    0.0, 0.0, 0.2, 1.0
), nrow = p, byrow = TRUE)
A2 <- matrix(c(
    1.8, 0.1, 0.0, 0.0,
    0.1, 0.7, 0.0, 0.0,
    0.0, 0.0, 0.8, -0.5,
    0.0, 0.0, -0.5, 1.5
), nrow = p, byrow = TRUE)
A3 <- matrix(c(
    0.9, -0.6, 0.2, 0.0,
    -0.6, 1.4, 0.0, 0.0,
    0.2, 0.0, 1.2, 0.6,
    0.0, 0.0, 0.6, 1.1
), nrow = p, byrow = TRUE)
X <- cbind(
    t(MASS::mvrnorm(60, mu = rep(0, p), Sigma = A1)),
    t(MASS::mvrnorm(60, mu = rep(0, p), Sigma = A2)),
    t(MASS::mvrnorm(60, mu = rep(0, p), Sigma = A3))
)
t0 <- proc.time()[["elapsed"]]
raw <- changepoints::BS.cov(X, 1L, n)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
elapsed <- proc.time()[["elapsed"]] - t0
if (is.null(trimmed$cpt_hat)) {
    cps_hat <- integer()
} else {
    cps_hat <- sort(as.integer(trimmed$cpt_hat[, 1]))
}

cat("n =", n, "\n")
cat("p =", p, "\n")
cat("tau =", tau, "\n")
cat("estimated changepoints =", paste(cps_hat, collapse = " "), "\n")
cat("raw nodes =", length(raw$S), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

set.seed(123)
n <- 120L
p <- 3L
lambda <- 0.5
zeta <- 20L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta1 <- c(1.5, 0.0, -1.0)
beta2 <- c(0.2, 1.8, -0.2)
beta3 <- c(-1.2, 0.5, 1.1)
y <- numeric(n)
y[1:40] <- X[1:40, ] %*% beta1 + rnorm(40, sd = 0.5)
y[41:80] <- X[41:80, ] %*% beta2 + rnorm(40, sd = 0.5)
y[81:120] <- X[81:120, ] %*% beta3 + rnorm(40, sd = 0.5)

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::DPDU.regression(y = y, X = X, lambda = lambda, zeta = zeta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("p =", p, "\n")
cat("lambda =", lambda, "\n")
cat("zeta =", zeta, "\n")
cat("estimated changepoints =", paste(as.numeric(ans$cpt), collapse = " "), "\n")
cat(sprintf("beta checksum = %.12f\n", sum(ans$beta_mat)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

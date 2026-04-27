#!/usr/bin/env Rscript

set.seed(123)
n <- 120L
p <- 3L
gamma_set <- c(1, 2, 3)
lambda_set <- c(0.2, 0.5, 1.0)
delta <- 5L
eps <- 0.001
true_cps <- c(40L, 80L)

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
beta1 <- c(1.5, 0.0, -1.0)
beta2 <- c(0.2, 1.8, -0.2)
beta3 <- c(-1.2, 0.5, 1.1)
y <- numeric(n)
y[1:40] <- X[1:40, , drop = FALSE] %*% beta1 + rnorm(40, sd = 0.5)
y[41:80] <- X[41:80, , drop = FALSE] %*% beta2 + rnorm(40, sd = 0.5)
y[81:120] <- X[81:120, , drop = FALSE] %*% beta3 + rnorm(40, sd = 0.5)

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::CV.search.DP.regression(y = y, X = X, gamma_set = gamma_set, lambda_set = lambda_set, delta = delta, eps = eps)
elapsed <- proc.time()[["elapsed"]] - t0
min_idx <- as.vector(arrayInd(which.min(ans$test_error), dim(ans$test_error)))

cat("n =", n, "\n")
cat("p =", p, "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("lambda_set =", paste(lambda_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
cat(sprintf("eps = %.3f\n", eps))
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
for (i in seq_along(lambda_set)) {
  for (j in seq_along(gamma_set)) {
    cat("gamma =", gamma_set[[j]], "\n")
    cat("lambda =", lambda_set[[i]], "\n")
    cat("estimated changepoints =", paste(unlist(ans$cpt_hat[j, i]), collapse = " "), "\n")
    cat("n changepoints =", unlist(ans$K_hat[j, i]), "\n")
    cat(sprintf("test error = %.12f\n", unlist(ans$test_error[j, i])))
    cat(sprintf("train error = %.12f\n", unlist(ans$train_error[j, i])))
  }
}
cat("best gamma =", gamma_set[[min_idx[1]]], "\n")
cat("best lambda =", lambda_set[[min_idx[2]]], "\n")
cat("best index =", paste(min_idx, collapse = " "), "\n")
cat("best changepoints =", paste(unlist(ans$cpt_hat[min_idx[1], min_idx[2]]), collapse = " "), "\n")
cat(sprintf("best test error = %.12f\n", unlist(ans$test_error[min_idx[1], min_idx[2]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

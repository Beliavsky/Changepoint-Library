#!/usr/bin/env Rscript

set.seed(0)
n <- 150L
b_value <- 3.5
true_cps <- c(100L)
y <- rnorm(n)
y[101:n] <- y[101:n] + 1.5
b_vec <- rep(b_value, n - 1L)

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::online.univar(y_vec = y, b_vec = b_vec)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("b_value =", b_value, "\n")
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
cat("estimated changepoint =", ans$cpt_hat, "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

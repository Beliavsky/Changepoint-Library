#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
set.seed(1)
n <- 200L
true_cp <- 100L
x <- c(rnorm(true_cp, 0, 1), rnorm(n - true_cp, 5, 1))

t0 <- proc.time()[["elapsed"]]
ans <- changepoint::cpt.mean(x, penalty = "SIC", method = "AMOC", class = FALSE, minseglen = 1L)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("penalty =", "SIC", "\n")
cat("minseglen =", 1L, "\n")
cat("true changepoint =", true_cp, "\n")
cat("estimated changepoint =", as.integer(ans[["cpt"]]), "\n")
cat(sprintf("conf.value = %.7f\n", as.numeric(ans[["conf.value"]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

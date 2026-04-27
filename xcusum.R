#!/usr/bin/env Rscript

set.seed(1)
n <- 200L
true_cp <- 100L
penalty_name <- "Manual"
penalty_value <- 0.5
minseglen <- 1L
x <- c(rnorm(true_cp, 0, 1), rnorm(n - true_cp, 5, 1))

t0 <- proc.time()[["elapsed"]]
est <- changepoint::cpt.mean(
    x,
    penalty = penalty_name,
    pen.value = penalty_value,
    method = "AMOC",
    test.stat = "CUSUM",
    class = FALSE,
    minseglen = minseglen
)
elapsed <- proc.time()[["elapsed"]] - t0
tmp <- changepoint:::single.mean.cusum.calc(x, extrainf = TRUE, minseglen = minseglen)
pen <- changepoint:::penalty_decision(
    penalty_name,
    penalty_value,
    n,
    diffparam = 1,
    asymcheck = "mean.cusum",
    method = "AMOC"
)

cat("n =", n, "\n")
cat("penalty =", penalty_name, "\n")
cat("pen.value =", penalty_value, "\n")
cat("minseglen =", minseglen, "\n")
cat("true changepoint =", true_cp, "\n")
cat("estimated changepoint =", as.integer(est), "\n")
cat(sprintf("test statistic = %.12f\n", as.numeric(tmp[["test statistic"]])))
cat(sprintf("effective penalty = %.12f\n", as.numeric(pen)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

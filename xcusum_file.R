#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcusum_data.txt"
penalty_name <- "Manual"
penalty_value <- 0.5
minseglen <- 1L

t0 <- proc.time()[["elapsed"]]
x <- read_series(data_file)
n <- length(x)
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

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("penalty =", penalty_name, "\n")
cat("pen.value =", penalty_value, "\n")
cat("minseglen =", minseglen, "\n")
cat("estimated changepoint =", as.integer(est), "\n")
cat(sprintf("test statistic = %.12f\n", as.numeric(tmp[["test statistic"]])))
cat(sprintf("effective penalty = %.12f\n", as.numeric(pen)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

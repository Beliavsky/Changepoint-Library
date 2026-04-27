#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xgamma_data.txt"
shape <- 2
penalty_name <- "MBIC"
minseglen <- 2L

t0 <- proc.time()[["elapsed"]]
x <- read_series(data_file)
n <- length(x)
est <- changepoint::cpt.meanvar(
    x,
    shape = shape,
    penalty = penalty_name,
    method = "AMOC",
    test.stat = "Gamma",
    class = FALSE,
    minseglen = minseglen
)
elapsed <- proc.time()[["elapsed"]] - t0
tmp <- changepoint:::single.meanvar.gamma.calc(x, shape = shape, extrainf = TRUE, minseglen = minseglen)
pen <- changepoint:::penalty_decision(
    penalty_name,
    0,
    n,
    diffparam = 1,
    asymcheck = "meanvar.gamma",
    method = "AMOC"
)
alt_use <- as.numeric(tmp[["alt"]]) + log(as.numeric(tmp[["cpt"]])) + log(n - as.numeric(tmp[["cpt"]]) + 1)

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("shape =", shape, "\n")
cat("penalty =", penalty_name, "\n")
cat("minseglen =", minseglen, "\n")
cat("estimated changepoint =", as.integer(est), "\n")
cat(sprintf("null = %.6f\n", as.numeric(tmp[["null"]])))
cat(sprintf("alt = %.6f\n", as.numeric(tmp[["alt"]])))
cat(sprintf("alt.mbic = %.6f\n", alt_use))
cat(sprintf("pen.value = %.6f\n", as.numeric(pen)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

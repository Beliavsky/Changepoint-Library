#!/usr/bin/env Rscript

set.seed(1)
n <- 200L
true_cp <- 100L
penalty_name <- "MBIC"
minseglen <- 1L
x <- c(rexp(true_cp, rate = 1 / 5), rexp(n - true_cp, rate = 1 / 12))

t0 <- proc.time()[["elapsed"]]
est <- changepoint::cpt.meanvar(
    x,
    penalty = penalty_name,
    method = "AMOC",
    test.stat = "Exponential",
    class = FALSE,
    minseglen = minseglen
)
elapsed <- proc.time()[["elapsed"]] - t0
tmp <- changepoint:::single.meanvar.exp.calc(x, extrainf = TRUE, minseglen = minseglen)
alt_use <- as.numeric(tmp[["alt"]]) + log(as.numeric(tmp[["cpt"]])) + log(n - as.numeric(tmp[["cpt"]]) + 1)

cat("n =", n, "\n")
cat("penalty =", penalty_name, "\n")
cat("minseglen =", minseglen, "\n")
cat("true changepoint =", true_cp, "\n")
cat("estimated changepoint =", as.integer(est[["cpt"]]), "\n")
cat(sprintf("p.value = %.7f\n", as.numeric(est[["p value"]])))
cat(sprintf("null = %.6f\n", as.numeric(tmp[["null"]])))
cat(sprintf("alt = %.6f\n", as.numeric(tmp[["alt"]])))
cat(sprintf("alt.mbic = %.6f\n", alt_use))
cat(sprintf("decision pen.value = %.6f\n", 0.0))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

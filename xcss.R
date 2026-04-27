#!/usr/bin/env Rscript

set.seed(1)
n <- 200L
true_cp <- 100L
x <- c(rnorm(true_cp, 0, 1), rnorm(n - true_cp, 0, 10))

t0 <- proc.time()[["elapsed"]]
ans <- changepoint::cpt.var(
  x,
  penalty = "Manual",
  pen.value = "log(2*log(n))",
  method = "AMOC",
  test.stat = "CSS",
  class = FALSE,
  minseglen = 2L
)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("penalty =", "Manual", "\n")
cat("pen.value =", "log(2*log(n))", "\n")
cat("minseglen =", 2L, "\n")
cat("true changepoint =", true_cp, "\n")
cat("estimated changepoint =", as.integer(ans[["cpt"]]), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

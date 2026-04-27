#!/usr/bin/env Rscript

set.seed(1)
n <- 200L
true_cp <- 100L
x <- 1:n
beta0 <- ifelse(x <= true_cp, 0, 50)
beta1 <- ifelse(x <= true_cp, 1, 0.25)
y <- beta0 + beta1 * x + rnorm(n)
data <- cbind(y, 1, x)

t0 <- proc.time()[["elapsed"]]
ans <- changepoint::cpt.reg(
  data,
  penalty = "SIC",
  method = "AMOC",
  test.stat = "Normal",
  class = FALSE,
  minseglen = 5L
)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("penalty =", "SIC", "\n")
cat("minseglen =", 5L, "\n")
cat("true changepoint =", true_cp, "\n")
cat("estimated changepoint =", as.integer(ans$cpts[["cpt"]]), "\n")
cat(sprintf("pen.value = %.5f\n", as.numeric(ans$pen.value)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

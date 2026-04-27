#!/usr/bin/env Rscript

n <- 120L
x <- seq(-1.0, 1.0, length.out = n)
set.seed(1127)
beta0 <- c(rep(-0.5, 40), rep(1.0, 40), rep(-1.2, 40))
beta1 <- c(rep(1.2, 40), rep(-0.8, 40), rep(0.9, 40))
y <- beta0 + beta1 * x + rnorm(n, sd = 0.18)
dat <- data.frame(y = y, x = x)

ans <- strucchange::efp(y ~ x, type = "OLS-CUSUM", data = dat)
proc <- as.numeric(ans$process)
cat("argmax abs process =", which.max(abs(proc)) - 1L, "\n")
cat(sprintf("max abs process = %.12f\n", max(abs(proc))))

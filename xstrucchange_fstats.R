#!/usr/bin/env Rscript

from <- 0.2
to <- 0.8
n <- 120L
x <- seq(-1.0, 1.0, length.out = n)
set.seed(913)
beta0 <- c(rep(-0.5, 40), rep(1.0, 40), rep(-1.2, 40))
beta1 <- c(rep(1.2, 40), rep(-0.8, 40), rep(0.9, 40))
y <- beta0 + beta1 * x + rnorm(n, sd = 0.18)
dat <- data.frame(y = y, x = x)

ans <- strucchange::Fstats(y ~ x, from = from, to = to, data = dat)
cat("breakpoint =", ans$breakpoint, "\n")
cat(sprintf("max F = %.12f\n", max(as.numeric(ans$Fstats))))

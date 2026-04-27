#!/usr/bin/env Rscript

h <- 20L
n <- 120L
x <- seq(-1.0, 1.0, length.out = n)
set.seed(641)
beta0 <- c(rep(-0.5, 40), rep(1.0, 40), rep(-1.2, 40))
beta1 <- c(rep(1.2, 40), rep(-0.8, 40), rep(0.9, 40))
y <- beta0 + beta1 * x + rnorm(n, sd = 0.18)
dat <- data.frame(y = y, x = x)

ans <- strucchange::breakpoints(y ~ x, h = h, data = dat)
sumry <- summary(ans)

cat("n =", n, "\n")
cat("p =", 1L, "\n")
cat("h =", h, "\n")
cat("best m =", which.min(as.numeric(sumry$RSS["BIC", ])) - 1L, "\n")
cat("estimated changepoints =", paste(ans$breakpoints, collapse = " "), "\n")

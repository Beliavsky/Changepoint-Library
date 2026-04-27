#!/usr/bin/env Rscript

set.seed(0)
n <- 1000L
h <- 120L
block_num <- 1L
guess_true <- 0.05
true_cps <- c(500L)
y <- c(rnorm(500, 0, 1), rnorm(500, 1, 1))
outlier_idx <- c(80L, 150L, 220L, 310L, 420L, 560L, 640L, 730L, 820L, 910L)
outlier_shift <- c(12, -11, 10, -13, 11, 12, -10, 13, -12, 11)
y[outlier_idx] <- y[outlier_idx] + outlier_shift
t_dat <- y[1:200]

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::aARC(y, t_dat = t_dat, guess_true = guess_true, h = h, block_num = block_num)
elapsed <- proc.time()[["elapsed"]] - t0

cat("n =", n, "\n")
cat("h =", h, "\n")
cat("block_num =", block_num, "\n")
cat("guess_true =", guess_true, "\n")
cat("true changepoints =", paste(true_cps, collapse = " "), "\n")
cat("estimated changepoints =", paste(ans, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

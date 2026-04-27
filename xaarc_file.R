#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xaarc_data.txt"
h <- 120L
block_num <- 1L
guess_true <- 0.05

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
t_dat <- y[1:200]
ans <- changepoints::aARC(y, t_dat = t_dat, guess_true = guess_true, h = h, block_num = block_num)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("t_dat_n =", length(t_dat), "\n")
cat("h =", h, "\n")
cat("block_num =", block_num, "\n")
cat("guess_true =", guess_true, "\n")
cat("estimated changepoints =", paste(ans, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

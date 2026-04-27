#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xarc_data.txt"
h <- 120L
block_num <- 1L
epsilon <- 0.1
gaussian <- TRUE

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
ans <- changepoints::ARC(y, h = h, block_num = block_num, epsilon = epsilon, gaussian = gaussian)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("h =", h, "\n")
cat("block_num =", block_num, "\n")
cat("epsilon =", epsilon, "\n")
cat("gaussian =", gaussian, "\n")
cat("estimated changepoints =", paste(ans, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

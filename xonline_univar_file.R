#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xonline_univar_data.txt"
threshold_file <- sub("\\.txt$", "_thresholds.txt", data_file)

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
b_vec <- read_series(threshold_file)
ans <- changepoints::online.univar(y_vec = y, b_vec = b_vec)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("threshold_file =", threshold_file, "\n")
cat("n =", length(y), "\n")
cat("estimated changepoint =", ans$cpt_hat, "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

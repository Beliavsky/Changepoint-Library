#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdp_univar_data.txt"
gamma <- 5
delta <- 5L

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
ans <- changepoints::DP.univar(y, gamma = gamma, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("gamma =", gamma, "\n")
cat("delta =", delta, "\n")
cat("estimated changepoints =", paste(ans$cpt, collapse = " "), "\n")
cat("n changepoints =", length(ans$cpt), "\n")
cat("partition last =", tail(ans$partition, 1), "\n")
cat(sprintf("yhat checksum = %.12f\n", sum(ans$yhat)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xamoc_data.txt"

t0 <- proc.time()[["elapsed"]]
x <- read_series(data_file)
n <- length(x)
ans <- changepoint::cpt.mean(x, penalty = "SIC", method = "AMOC", class = FALSE, minseglen = 1L)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("penalty =", "SIC", "\n")
cat("minseglen =", 1L, "\n")
cat("estimated changepoint =", as.integer(ans[["cpt"]]), "\n")
cat(sprintf("conf.value = %.7f\n", as.numeric(ans[["conf.value"]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

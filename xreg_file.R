#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xreg_data.txt"

t0 <- proc.time()[["elapsed"]]
data <- read_matrix(data_file)
ans <- changepoint::cpt.reg(
  data,
  penalty = "SIC",
  method = "AMOC",
  test.stat = "Normal",
  class = FALSE,
  minseglen = 5L
)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", nrow(data), "\n")
cat("penalty =", "SIC", "\n")
cat("minseglen =", 5L, "\n")
cat("estimated changepoint =", as.integer(ans$cpts[["cpt"]]), "\n")
cat(sprintf("pen.value = %.5f\n", as.numeric(ans$pen.value)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

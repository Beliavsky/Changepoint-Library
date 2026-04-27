#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcss_data.txt"

t0 <- proc.time()[["elapsed"]]
x <- read_series(data_file)
ans <- changepoint::cpt.var(
  x,
  penalty = "Manual",
  pen.value = "log(2*log(n))",
  method = "AMOC",
  test.stat = "CSS",
  class = FALSE,
  minseglen = 2L
)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(x), "\n")
cat("penalty =", "Manual", "\n")
cat("pen.value =", "log(2*log(n))", "\n")
cat("minseglen =", 2L, "\n")
cat("estimated changepoint =", as.integer(ans[["cpt"]]), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

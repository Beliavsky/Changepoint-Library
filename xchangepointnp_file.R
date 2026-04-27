#!/usr/bin/env Rscript

library(changepoint.np)

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xchangepointnp_data.txt"
penalty <- "MBIC"
minseglen <- 2L
nquantiles <- 10L

x <- scan(data_file, quiet = TRUE, comment.char = "#")
pen.value <- changepoint:::penalty_decision(penalty, 0, length(x), 1, method = "PELT")
data_input_np <- get("data_input", asNamespace("changepoint.np"))

t0 <- proc.time()[["elapsed"]]
out <- data_input_np(
  data = x,
  method = "PELT",
  pen.value = pen.value,
  costfunc = "nonparametric.ed.mbic",
  minseglen = minseglen,
  nquantiles = nquantiles
)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(x), "\n")
cat("penalty =", penalty, "\n")
cat(sprintf("pen.value = %.12f\n", pen.value))
cat("minseglen =", minseglen, "\n")
cat("nquantiles =", nquantiles, "\n")
cat("cpts =", paste(out[[2]], collapse = " "), "\n")
cat(sprintf("lastchangecpts checksum = %.12f\n", sum(out[[1]])))
cat(sprintf("lastchangelike checksum = %.12f\n", sum(out[[3]])))
cat(sprintf("numchangecpts checksum = %.12f\n", sum(out[[4]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

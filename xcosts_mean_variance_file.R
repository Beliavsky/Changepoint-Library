#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcosts_mean_variance_data.txt"
n_bkps <- 3L
min_size <- 30L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
models <- list(
  list(model = "l1", note = "robust location shifts"),
  list(model = "l2", note = "mean shifts"),
  list(model = "normal", note = "mean and variance shifts"),
  list(model = "rbf", note = "general distribution shifts")
)

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("n_bkps =", n_bkps, "\n")
cat("min_size =", min_size, "\n")
cat("\n")

for (item in models) {
  t1 <- proc.time()[["elapsed"]]
  bkps <- solve_dynp_cost_1d(signal, n_bkps, min_size, item$model)
  elapsed <- proc.time()[["elapsed"]] - t1
  cat("model =", item$model, "\n")
  cat("note =", item$note, "\n")
  cat("estimated bkps =", paste(bkps, collapse = " "), "\n")
  cat(sprintf("elapsed seconds = %.3f\n", elapsed))
  cat("\n")
}

cat(sprintf("wall time elapsed (s) = %.3f\n", proc.time()[["elapsed"]] - t0))

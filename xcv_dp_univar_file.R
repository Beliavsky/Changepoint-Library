#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcv_dp_univar_data.txt"
gamma_set <- 1:8
delta <- 5L

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
ans <- changepoints::CV.search.DP.univar(y, gamma_set = gamma_set, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0
min_idx <- which.min(unlist(ans$test_error))

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
for (j in seq_along(gamma_set)) {
  cat("gamma =", gamma_set[[j]], "\n")
  cat("estimated changepoints =", paste(ans$cpt_hat[[j]], collapse = " "), "\n")
  cat("n changepoints =", ans$K_hat[[j]], "\n")
  cat(sprintf("test error = %.12f\n", ans$test_error[[j]]))
  cat(sprintf("train error = %.12f\n", ans$train_error[[j]]))
}
cat("best gamma =", gamma_set[[min_idx]], "\n")
cat("best index =", min_idx, "\n")
cat("best changepoints =", paste(ans$cpt_hat[[min_idx]], collapse = " "), "\n")
cat(sprintf("best test error = %.12f\n", ans$test_error[[min_idx]]))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

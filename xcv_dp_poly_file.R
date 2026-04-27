#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcv_dp_poly_data.txt"
r <- 2L
gamma_set <- c(0.5, 1, 2, 3)
delta <- 5L

y <- as.numeric(read.table(data_file, comment.char = "#")[, 1])

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::CV.search.DP.poly(y = y, r = r, gamma_set = gamma_set, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0
min_idx <- which.min(unlist(ans$test_error))

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("r =", r, "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
for (j in seq_along(gamma_set)) {
  cat("gamma =", gamma_set[[j]], "\n")
  cat("estimated changepoints =", paste(unlist(ans$cpt_hat[j]), collapse = " "), "\n")
  cat("n changepoints =", unlist(ans$K_hat[j]), "\n")
  cat(sprintf("test error = %.12f\n", unlist(ans$test_error[j])))
  cat(sprintf("train error = %.12f\n", unlist(ans$train_error[j])))
}
cat("best gamma =", gamma_set[[min_idx]], "\n")
cat("best index =", min_idx, "\n")
cat("best changepoints =", paste(unlist(ans$cpt_hat[min_idx]), collapse = " "), "\n")
cat(sprintf("best test error = %.12f\n", unlist(ans$test_error[min_idx])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

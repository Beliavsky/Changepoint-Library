#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcv_dp_var1_data.txt"
gamma_set <- c(0.01, 0.05, 0.1)
lambda_set <- c(0, 0.05, 0.1)
delta <- 5L
eps <- 0.001

DATA <- t(as.matrix(read.table(data_file, comment.char = "#")))

t0 <- proc.time()[["elapsed"]]
ans <- changepoints:::CV.search.DP.VAR1(DATA, gamma_set = gamma_set, lambda_set = lambda_set, delta = delta, eps = eps)
elapsed <- proc.time()[["elapsed"]] - t0
min_idx <- as.vector(arrayInd(which.min(ans$test_error), dim(ans$test_error)))

cat("file =", data_file, "\n")
cat("p =", nrow(DATA), "\n")
cat("n =", ncol(DATA), "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("lambda_set =", paste(lambda_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
cat(sprintf("eps = %.3f\n", eps))
for (i in seq_along(lambda_set)) {
  for (j in seq_along(gamma_set)) {
    cat("gamma =", gamma_set[[j]], "\n")
    cat("lambda =", lambda_set[[i]], "\n")
    cat("estimated changepoints =", paste(unlist(ans$cpt_hat[j, i]), collapse = " "), "\n")
    cat("n changepoints =", unlist(ans$K_hat[j, i]), "\n")
    cat(sprintf("test error = %.12f\n", unlist(ans$test_error[j, i])))
    cat(sprintf("train error = %.12f\n", unlist(ans$train_error[j, i])))
  }
}
cat("best gamma =", gamma_set[[min_idx[1]]], "\n")
cat("best lambda =", lambda_set[[min_idx[2]]], "\n")
cat("best index =", paste(min_idx, collapse = " "), "\n")
cat("best changepoints =", paste(unlist(ans$cpt_hat[min_idx[1], min_idx[2]]), collapse = " "), "\n")
cat(sprintf("best test error = %.12f\n", unlist(ans$test_error[min_idx[1], min_idx[2]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

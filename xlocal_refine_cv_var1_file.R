#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcv_dp_var1_data.txt"
gamma_set <- c(0.01, 0.05, 0.1)
lambda_set <- c(0, 0.05, 0.1)
delta <- 5L
eps <- 0.001
zeta_set <- c(0.01, 0.05, 0.1, 0.2)
delta_local <- 5L

DATA <- t(as.matrix(read.table(data_file, comment.char = "#")))

t0 <- proc.time()[["elapsed"]]
ans_cv <- changepoints:::CV.search.DP.VAR1(DATA, gamma_set = gamma_set, lambda_set = lambda_set, delta = delta, eps = eps)
min_idx <- as.vector(arrayInd(which.min(ans_cv$test_error), dim(ans_cv$test_error)))
cpt_init <- unlist(ans_cv$cpt_hat[min_idx[1], min_idx[2]])
ans <- changepoints:::local.refine.CV.VAR1(cpt_init = cpt_init, DATA = DATA, zeta_set = zeta_set, delta_local = delta_local)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("p =", nrow(DATA), "\n")
cat("n =", ncol(DATA), "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("lambda_set =", paste(lambda_set, collapse = " "), "\n")
cat("zeta_set =", paste(zeta_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
cat("delta_local =", delta_local, "\n")
cat("best gamma =", gamma_set[[min_idx[1]]], "\n")
cat("best lambda =", lambda_set[[min_idx[2]]], "\n")
cat("initial changepoints =", paste(cpt_init, collapse = " "), "\n")
cat("refined changepoints =", paste(ans$cpt_hat, collapse = " "), "\n")
cat("zeta =", ans$zeta, "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

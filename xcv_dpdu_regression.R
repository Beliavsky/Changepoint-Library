#!/usr/bin/env Rscript

data_file <- "xcv_dpdu_regression_data.txt"
lambda_set <- c(0.2, 0.5, 1.0)
zeta_set <- c(10L, 20L, 30L)

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
X <- dat[, -1, drop = FALSE]

t0 <- proc.time()[["elapsed"]]
ans <- changepoints:::CV.search.DPDU.regression(y = y, X = X, lambda_set = lambda_set, zeta_set = zeta_set)
elapsed <- proc.time()[["elapsed"]] - t0
min_idx <- as.vector(arrayInd(which.min(ans$test_error), dim(ans$test_error)))

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("p =", ncol(X), "\n")
cat("lambda_set =", paste(lambda_set, collapse = " "), "\n")
cat("zeta_set =", paste(zeta_set, collapse = " "), "\n")
for (i in seq_along(lambda_set)) {
  for (j in seq_along(zeta_set)) {
    cat("lambda =", lambda_set[[i]], "\n")
    cat("zeta =", zeta_set[[j]], "\n")
    cat("estimated changepoints =", paste(unlist(ans$cpt_hat[j, i]), collapse = " "), "\n")
    cat("n changepoints =", unlist(ans$K_hat[j, i]), "\n")
    cat(sprintf("test error = %.12f\n", unlist(ans$test_error[j, i])))
    cat(sprintf("train error = %.12f\n", unlist(ans$train_error[j, i])))
  }
}
cat("best lambda =", lambda_set[[min_idx[2]]], "\n")
cat("best zeta =", zeta_set[[min_idx[1]]], "\n")
cat("best index =", paste(min_idx, collapse = " "), "\n")
cat("best changepoints =", paste(unlist(ans$cpt_hat[min_idx[1], min_idx[2]]), collapse = " "), "\n")
cat(sprintf("best test error = %.12f\n", unlist(ans$test_error[min_idx[1], min_idx[2]])))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

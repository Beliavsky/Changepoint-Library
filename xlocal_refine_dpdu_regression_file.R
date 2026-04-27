#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xdpdu_regression_data.txt"
lambda <- 0.5
zeta <- 20L
w <- 0.9

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
X <- dat[, -1, drop = FALSE]

t0 <- proc.time()[["elapsed"]]
init <- changepoints::DPDU.regression(y = y, X = X, lambda = lambda, zeta = zeta)
beta_hat <- init$beta_mat[, c(init$cpt, nrow(X)), drop = FALSE]
refined <- changepoints:::local.refine.DPDU.regression(init$cpt, beta_hat, y = y, X = X, w = w)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("p =", ncol(X), "\n")
cat("lambda =", lambda, "\n")
cat("zeta =", zeta, "\n")
cat("w =", w, "\n")
cat("initial changepoints =", paste(init$cpt, collapse = " "), "\n")
cat("refined changepoints =", paste(refined, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

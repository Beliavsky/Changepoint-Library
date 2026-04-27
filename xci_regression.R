#!/usr/bin/env Rscript

data_file <- "xdpdu_regression_data.txt"
u_file <- "xci_regression_u.txt"
lambda <- 0.5
zeta <- 20L
w <- 0.9
B <- 200L
alpha_vec <- c(0.05, 0.10)

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
X <- dat[, -1, drop = FALSE]
n <- length(y)
M <- n

t0 <- proc.time()[["elapsed"]]
init <- changepoints::DPDU.regression(y = y, X = X, lambda = lambda, zeta = zeta)
beta_hat <- init$beta_mat[, c(init$cpt, nrow(X)), drop = FALSE]
cpt_lr <- changepoints:::local.refine.DPDU.regression(init$cpt, beta_hat, y = y, X = X, w = w)
ci <- changepoints::CI.regression(init$cpt, cpt_lr, beta_hat, y, X, w = w, B = B, M = M, alpha_vec = alpha_vec)
interval_refine_mat <- changepoints:::trim_interval(n = n, init$cpt, w = w)
block_size <- ceiling((min(floor(interval_refine_mat[, 2]) - ceiling(interval_refine_mat[, 1])))^(2 / 5) / 2)
lrv_hat <- changepoints:::LRV.regression(init$cpt, beta_hat, y, X, w = w, block_size = block_size)
kappa2_hat <- apply(beta_hat[, -ncol(beta_hat), drop = FALSE] - beta_hat[, -1, drop = FALSE], MARGIN = 2, crossprod)
drift_hat <- diag(t(beta_hat[, -ncol(beta_hat), drop = FALSE] - beta_hat[, -1, drop = FALSE]) %*%
  (t(cbind(rep(1, n), X)) %*% cbind(rep(1, n), X)) %*%
  (beta_hat[, -ncol(beta_hat), drop = FALSE] - beta_hat[, -1, drop = FALSE]) / (n * kappa2_hat))
u_mat <- matrix(NA_real_, nrow = length(init$cpt), ncol = B)
for (i in seq_along(init$cpt)) {
  for (b in seq_len(B)) {
    set.seed(12345 + b + 10000 * i)
    u_mat[i, b] <- seq(-M, M)[which.min(changepoints:::simu.2BM_Drift(M, drift_hat[i], lrv_hat[i]))]
  }
}
elapsed <- proc.time()[["elapsed"]] - t0

write.table(u_mat, file = u_file, row.names = FALSE, col.names = FALSE)
cat("file =", data_file, "\n")
cat("u_file =", u_file, "\n")
cat("n =", length(y), "\n")
cat("p =", ncol(X), "\n")
cat("lambda =", lambda, "\n")
cat("zeta =", zeta, "\n")
cat("w =", w, "\n")
cat("B =", B, "\n")
cat("M =", M, "\n")
cat("alpha_vec =", paste(sprintf('%.2f', alpha_vec), collapse = " "), "\n")
cat("initial changepoints =", paste(init$cpt, collapse = " "), "\n")
cat("refined changepoints =", paste(cpt_lr, collapse = " "), "\n")
cat("block_size =", block_size, "\n")
cat("lrv_hat =", paste(sprintf('%.12f', lrv_hat), collapse = " "), "\n")
cat("kappa2_hat =", paste(sprintf('%.12f', kappa2_hat), collapse = " "), "\n")
cat("drift_hat =", paste(sprintf('%.12f', drift_hat), collapse = " "), "\n")
for (j in seq_along(alpha_vec)) {
  cat("alpha =", sprintf("%.2f", alpha_vec[j]), "\n")
  for (i in seq_along(init$cpt)) {
    cat("ci =", i, sprintf("%.0f", ci[i, 1, j]), sprintf("%.0f", ci[i, 2, j]), "\n")
  }
}
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

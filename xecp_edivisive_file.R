#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xecp_edivisive_data.txt"
k <- 2L
min_size <- 10L
alpha <- 1

X <- as.matrix(read.table(data_file, comment.char = "#"))
t0 <- proc.time()[["elapsed"]]
fit <- ecp::e.divisive(X, k = k, min.size = min_size, alpha = alpha)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", nrow(X), "\n")
cat("p =", ncol(X), "\n")
cat("k =", k, "\n")
cat("min.size =", min_size, "\n")
cat("alpha =", alpha, "\n")
cat("k.hat =", fit$k.hat, "\n")
cat("estimates =", paste(fit$estimates, collapse = " "), "\n")
cat("order.found =", paste(fit$order.found, collapse = " "), "\n")
cat(sprintf("cluster checksum = %.12f\n", sum(fit$cluster)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

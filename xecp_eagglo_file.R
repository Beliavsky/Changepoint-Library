#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xecp_eagglo_data.txt"
alpha <- 1

X <- as.matrix(read.table(data_file, comment.char = "#"))
t0 <- proc.time()[["elapsed"]]
fit <- ecp::e.agglo(X, alpha = alpha, penalty = function(cps) 0)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", nrow(X), "\n")
cat("p =", ncol(X), "\n")
cat("alpha =", alpha, "\n")
cat("estimates =", paste(fit$estimates, collapse = " "), "\n")
cat(sprintf("cluster checksum = %.12f\n", sum(fit$cluster)))
cat("fit length =", length(fit$fit), "\n")
cat(sprintf("fit checksum = %.12f\n", sum(fit$fit)))
cat("progression dim =", nrow(fit$progression), ncol(fit$progression), "\n")
cat(sprintf("progression checksum = %.12f\n", sum(fit$progression, na.rm = TRUE)))
cat(sprintf("merged checksum = %.12f\n", sum(fit$merged, na.rm = TRUE)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

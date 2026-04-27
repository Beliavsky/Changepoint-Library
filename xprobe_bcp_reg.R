#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bcp))

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xbcp_reg_data.txt"
burnin <- 50L
mcmc <- 200L
p0 <- 0.2
w0 <- c(0.2, 0.2)

dat <- as.matrix(read.table(data_file, comment.char = "#"))
y <- dat[, 1]
x <- dat[, 2]
RNGkind("Wichmann-Hill")
set.seed(123)

fit <- bcp(y, x, burnin = burnin, mcmc = mcmc, p0 = p0, w0 = w0, return.mcmc = TRUE)
rho_mat <- fit$mcmc.rhos[, (burnin + 1):(burnin + mcmc), drop = FALSE]
rho_counts <- rowSums(rho_mat)
pp <- fit$posterior.prob
keep <- !is.na(pp)
idx <- which(rho_counts[keep] != 0)

cat("probe pp checksum =", sprintf("%.15f", sum(pp[keep])), "\n")
cat("probe rho counts checksum =", sprintf("%.15f", sum(rho_counts)), "\n")
cat("probe top idx =", paste(order(pp[keep], decreasing = TRUE)[1:5], collapse = " "), "\n")
cat("probe blocks head =", paste(fit$blocks[1:20], collapse = " "), "\n")
cat("probe rho counts head =", paste(rho_counts[1:20], collapse = " "), "\n")
cat("probe rho counts around 40 =", paste(rho_counts[35:45], collapse = " "), "\n")
cat("probe rho counts around 80 =", paste(rho_counts[75:85], collapse = " "), "\n")
cat("probe rho nz =", paste(paste(idx, rho_counts[idx], sep = ":"), collapse = " "), "\n")

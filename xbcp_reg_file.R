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
seed_vals <- .Random.seed[2:4]

t0 <- proc.time()[["elapsed"]]
fit <- bcp(y, x, burnin = burnin, mcmc = mcmc, p0 = p0, w0 = w0)
elapsed <- proc.time()[["elapsed"]] - t0

pm <- fit$posterior.mean[, 1]
pp <- fit$posterior.prob
keep <- !is.na(pp)
ord <- order(pp[keep], decreasing = TRUE)[seq_len(5)]
top_idx <- which(keep)[ord]

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("burnin =", burnin, "\n")
cat("mcmc =", mcmc, "\n")
cat(sprintf("p0 = %.12f\n", p0))
cat("w0 =", paste(sprintf("%.12f", w0), collapse = " "), "\n")
cat("rng = Wichmann-Hill\n")
cat("seed =", paste(seed_vals, collapse = " "), "\n")
cat(sprintf("pm checksum = %.15f\n", sum(pm)))
cat(sprintf("pp checksum = %.15f\n", sum(pp[keep])))
cat("top idx =", paste(top_idx, collapse = " "), "\n")
cat("top val =", paste(sprintf("%.15f", pp[top_idx]), collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

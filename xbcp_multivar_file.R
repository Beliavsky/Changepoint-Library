#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bcp))

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xbcp_multivar_data.txt"
burnin <- 50L
mcmc <- 200L
p0 <- 0.2
w0 <- 0.2

y <- as.matrix(read.table(data_file, comment.char = "#"))
RNGkind("Wichmann-Hill")
set.seed(123)
seed_vals <- .Random.seed[2:4]

t0 <- proc.time()[["elapsed"]]
fit <- bcp(y, burnin = burnin, mcmc = mcmc, p0 = p0, w0 = w0)
elapsed <- proc.time()[["elapsed"]] - t0

pm <- fit$posterior.mean
pp <- fit$posterior.prob
keep <- !is.na(pp)
ord <- order(pp[keep], decreasing = TRUE)[seq_len(5)]
top_idx <- which(keep)[ord]

cat("file =", data_file, "\n")
cat("n =", nrow(y), "\n")
cat("p =", ncol(y), "\n")
cat("burnin =", burnin, "\n")
cat("mcmc =", mcmc, "\n")
cat(sprintf("p0 = %.12f\n", p0))
cat(sprintf("w0 = %.12f\n", w0))
cat("rng = Wichmann-Hill\n")
cat("seed =", paste(seed_vals, collapse = " "), "\n")
cat(sprintf("pm checksum = %.15f\n", sum(pm)))
cat(sprintf("pp checksum = %.15f\n", sum(pp[keep])))
cat("top idx =", paste(top_idx, collapse = " "), "\n")
cat("top val =", paste(sprintf("%.15f", pp[top_idx]), collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

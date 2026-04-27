#!/usr/bin/env Rscript

data_file <- "xdp_poly_data.txt"
r <- 2L
gamma <- 1
delta <- 5L

y <- as.numeric(read.table(data_file, comment.char = "#")[, 1])

t0 <- proc.time()[["elapsed"]]
ans <- changepoints::DP.poly(y = y, r = r, gamma = gamma, delta = delta)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("r =", r, "\n")
cat("gamma =", gamma, "\n")
cat("delta =", delta, "\n")
cat("estimated changepoints =", paste(ans$cpt, collapse = " "), "\n")
cat(sprintf("yhat checksum = %.12f\n", sum(ans$yhat)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

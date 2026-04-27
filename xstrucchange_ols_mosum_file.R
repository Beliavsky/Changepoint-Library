#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_ols_mosum_data.txt"
h <- 0.25

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")

t0 <- proc.time()[["elapsed"]]
ans <- strucchange::efp(y ~ x, type = "OLS-MOSUM", h = h, data = dat)
elapsed <- proc.time()[["elapsed"]] - t0

proc <- as.numeric(ans$process)
imax <- which.max(abs(proc))
cat("file =", data_file, "\n")
cat("n =", nrow(dat), "\n")
cat("p =", 1L, "\n")
cat("type =", "OLS-MOSUM", "\n")
cat("h =", h, "\n")
cat("argmax abs process =", imax - 1L, "\n")
cat(sprintf("max abs process = %.12f\n", max(abs(proc))))
cat(sprintf("process checksum = %.12f\n", sum(proc)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

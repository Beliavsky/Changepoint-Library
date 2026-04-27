#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_mefp_data.txt"
history_n <- 50L
alpha <- 0.05

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")
hist_dat <- dat[1:history_n, , drop = FALSE]

t0 <- proc.time()[["elapsed"]]
me0 <- strucchange::mefp(y ~ x, type = "OLS-CUSUM", data = hist_dat, alpha = alpha)
ans <- strucchange::monitor(me0, data = dat, verbose = FALSE)
elapsed <- proc.time()[["elapsed"]] - t0

proc <- as.numeric(ans$process)
bd <- as.numeric(ans$border((history_n + 1L):nrow(dat)))
cat("file =", data_file, "\n")
cat("n =", nrow(dat), "\n")
cat("history n =", history_n, "\n")
cat("p =", 1L, "\n")
cat("type =", "mefp OLS-CUSUM", "\n")
cat("alpha =", alpha, "\n")
cat(sprintf("critval = %.12f\n", ans$critval))
cat("breakpoint =", ans$breakpoint, "\n")
cat("argmax abs process =", which.max(abs(proc)) - 1L, "\n")
cat(sprintf("max abs process = %.12f\n", max(abs(proc))))
cat(sprintf("process checksum = %.12f\n", sum(proc)))
cat(sprintf("boundary checksum = %.12f\n", sum(bd)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

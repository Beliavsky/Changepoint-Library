#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_me_data.txt"
h <- 0.25

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")

t0 <- proc.time()[["elapsed"]]
ans <- strucchange::efp(y ~ x, type = "ME", h = h, data = dat)
elapsed <- proc.time()[["elapsed"]] - t0

proc <- ans[["process"]]
idx <- which(abs(proc) == max(abs(proc)), arr.ind = TRUE)[1, ]
cat("file =", data_file, "\n")
cat("n =", nrow(dat), "\n")
cat("p =", 1L, "\n")
cat("type =", "ME", "\n")
cat("h =", h, "\n")
cat("nrow process =", nrow(proc), "\n")
cat("ncol process =", ncol(proc), "\n")
cat("imax row =", idx[1], "\n")
cat("imax col =", idx[2], "\n")
cat(sprintf("max abs process = %.12f\n", max(abs(proc))))
cat(sprintf("process checksum = %.12f\n", sum(proc)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

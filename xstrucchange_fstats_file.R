#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_fstats_data.txt"
from <- 0.2
to <- 0.8

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")

t0 <- proc.time()[["elapsed"]]
ans <- strucchange::Fstats(y ~ x, from = from, to = to, data = dat)
elapsed <- proc.time()[["elapsed"]] - t0

fst <- as.numeric(ans$Fstats)
cat("file =", data_file, "\n")
cat("n =", nrow(dat), "\n")
cat("p =", 1L, "\n")
cat("from =", from, "\n")
cat("to =", to, "\n")
cat("breakpoint =", ans$breakpoint, "\n")
cat(sprintf("max F = %.12f\n", max(fst)))
cat(sprintf("min RSS = %.12f\n", ans$RSS))
cat(sprintf("Fstats checksum = %.12f\n", sum(fst)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

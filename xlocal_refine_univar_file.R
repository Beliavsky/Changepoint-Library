#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xwbs_data.txt"
interval_file <- sub("\\.txt$", "_intervals.txt", data_file)
delta <- 5L
tau <- 4

t0 <- proc.time()[["elapsed"]]
y <- read_series(data_file)
intervals <- read.table(interval_file, comment.char = "#")
raw <- changepoints::WBS.univar(y, 1L, length(y), intervals[[1]], intervals[[2]], delta = delta)
trimmed <- changepoints::thresholdBS(raw, tau = tau)
if (is.null(trimmed$cpt_hat)) {
    cps_init <- integer()
    cps_refined <- integer()
} else {
    cps_init <- sort(as.integer(trimmed$cpt_hat[, 1]))
    cps_refined <- changepoints::local.refine.univar(cps_init, y)
}
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("interval_file =", interval_file, "\n")
cat("n =", length(y), "\n")
cat("M =", nrow(intervals), "\n")
cat("delta =", delta, "\n")
cat("tau =", tau, "\n")
cat("initial changepoints =", paste(cps_init, collapse = " "), "\n")
cat("refined changepoints =", paste(cps_refined, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

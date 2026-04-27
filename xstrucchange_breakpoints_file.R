#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_breakpoints_data.txt"
h <- 20L

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")

t0 <- proc.time()[["elapsed"]]
ans <- strucchange::breakpoints(y ~ x, h = h, data = dat)
sumry <- summary(ans)
elapsed <- proc.time()[["elapsed"]] - t0

bic_vals <- as.numeric(sumry$RSS["BIC", ])
rss_vals <- as.numeric(sumry$RSS["RSS", ])
best_m <- which.min(bic_vals) - 1L
cpt <- ans$breakpoints
if (all(is.na(cpt))) cpt <- integer(0)

cat("file =", data_file, "\n")
cat("n =", nrow(dat), "\n")
cat("p =", 1L, "\n")
cat("h =", h, "\n")
cat("max breaks =", length(bic_vals) - 1L, "\n")
cat("best m =", best_m, "\n")
cat("estimated changepoints =", paste(cpt, collapse = " "), "\n")
cat("RSS =", paste(sprintf("%.12f", rss_vals), collapse = " "), "\n")
cat("BIC =", paste(sprintf("%.12f", bic_vals), collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

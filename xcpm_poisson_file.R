#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcpm_poisson_data.txt"
threshold_file <- "xcpm_poisson_thresholds.txt"
batch_threshold_file <- "xcpm_poisson_batch_threshold.txt"
ARL0 <- NA
startup <- 20
alpha <- NA

x <- scan(data_file, quiet = TRUE, comment.char = "#")
thresholds <- rep(999999, length(x))
batch_threshold <- -1
write.table(thresholds, file = threshold_file, row.names = FALSE, col.names = FALSE)
write.table(batch_threshold, file = batch_threshold_file, row.names = FALSE, col.names = FALSE)

t0 <- proc.time()[["elapsed"]]
det <- cpm::detectChangePoint(x, "Poisson", ARL0 = ARL0, startup = startup)
proc <- cpm::processStream(x, "Poisson", ARL0 = ARL0, startup = startup)
bat <- cpm::detectChangePointBatch(x, "Poisson", alpha = alpha)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(x), "\n")
cat("type = Poisson\n")
cat("ARL0 = NA\n")
cat("startup =", startup, "\n")
cat("alpha = NA\n")
cat("detect changePoint =", det$changePoint, "\n")
cat("detect detectionTime =", det$detectionTime, "\n")
cat(sprintf("detect Ds checksum = %.12f\n", sum(det$Ds)))
cat("process changePoints =", paste(proc$changePoints, collapse = " "), "\n")
cat("process detectionTimes =", paste(proc$detectionTimes, collapse = " "), "\n")
cat("batch changePoint =", bat$changePoint, "\n")
cat(sprintf("batch threshold = %.12f\n", batch_threshold))
cat(sprintf("batch Ds checksum = %.12f\n", sum(bat$Ds)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcpm_fet_data.txt"
threshold_file <- "xcpm_fet_thresholds.txt"
batch_threshold_file <- "xcpm_fet_batch_threshold.txt"
ARL0 <- 500
startup <- 20
alpha <- 0.05
lambda <- 0.3

x <- scan(data_file, quiet = TRUE, comment.char = "#")
thresholds <- cpm:::loadThresholds("FET", ARL0, length(x), lambda = lambda)
batch_threshold <- cpm::getBatchThreshold("FET", alpha, length(x), lambda = lambda)
write.table(thresholds, file = threshold_file, row.names = FALSE, col.names = FALSE)
write.table(batch_threshold, file = batch_threshold_file, row.names = FALSE, col.names = FALSE)

t0 <- proc.time()[["elapsed"]]
det <- cpm::detectChangePoint(x, "FET", ARL0 = ARL0, startup = startup, lambda = lambda)
proc <- cpm::processStream(x, "FET", ARL0 = ARL0, startup = startup, lambda = lambda)
bat <- cpm::detectChangePointBatch(x, "FET", alpha = alpha, lambda = lambda)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(x), "\n")
cat("type = FET\n")
cat("ARL0 =", ARL0, "\n")
cat("startup =", startup, "\n")
cat("alpha =", alpha, "\n")
cat("lambda =", lambda, "\n")
cat("detect changePoint =", det$changePoint, "\n")
cat("detect detectionTime =", det$detectionTime, "\n")
cat(sprintf("detect Ds checksum = %.12f\n", sum(det$Ds)))
cat("process changePoints =", paste(proc$changePoints, collapse = " "), "\n")
cat("process detectionTimes =", paste(proc$detectionTimes, collapse = " "), "\n")
cat("batch changePoint =", bat$changePoint, "\n")
cat(sprintf("batch threshold = %.12f\n", bat$threshold))
cat(sprintf("batch Ds checksum = %.12f\n", sum(bat$Ds)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

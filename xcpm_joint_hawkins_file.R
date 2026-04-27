#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xcpm_joint_hawkins_data.txt"
threshold_file <- "xcpm_joint_hawkins_thresholds.txt"
batch_threshold_file <- "xcpm_joint_hawkins_batch_threshold.txt"
ARL0 <- 500
startup <- 20
alpha <- NA

x <- scan(data_file, quiet = TRUE, comment.char = "#")
thresholds <- cpm:::loadThresholds("JointHawkins", ARL0, length(x), lambda = NA)
batch_threshold <- -1
write.table(thresholds, file = threshold_file, row.names = FALSE, col.names = FALSE)
write.table(batch_threshold, file = batch_threshold_file, row.names = FALSE, col.names = FALSE)

t0 <- proc.time()[["elapsed"]]
det <- cpm::detectChangePoint(x, "JointHawkins", ARL0 = ARL0, startup = startup)
bat <- cpm::detectChangePointBatch(x, "JointHawkins", alpha = alpha)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(x), "\n")
cat("type = JointHawkins\n")
cat("ARL0 =", ARL0, "\n")
cat("startup =", startup, "\n")
cat("alpha = NA\n")
cat("detect changePoint =", det$changePoint, "\n")
cat("detect detectionTime =", det$detectionTime, "\n")
cat(sprintf("detect Ds checksum = %.12f\n", sum(det$Ds)))
cat("batch changePoint =", bat$changePoint, "\n")
cat("batch threshold = -1.000000000000\n")
cat(sprintf("batch Ds checksum = %.12f\n", sum(bat$Ds)))
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

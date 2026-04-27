#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmcp_demo_data.txt"
draws_file <- "xmcp_demo_draws.txt"

e <- new.env()
lazyLoad(file.path(find.package("mcp"), "data", "Rdata"), e)
fit <- get("demo_fit", envir = e)
dat <- fit[["data"]]
draws <- do.call(rbind, lapply(fit[["mcmc_post"]], as.matrix))

write.table(dat, file = data_file, row.names = FALSE, col.names = FALSE)
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

time <- dat[["time"]]
cp1 <- draws[, "cp_1"]
cp2 <- draws[, "cp_2"]
int1 <- draws[, "int_1"]
int3 <- draws[, "int_3"]
time2 <- draws[, "time_2"]
time3 <- draws[, "time_3"]
sigma1 <- draws[, "sigma_1"]

mu <- sapply(time, function(t) {
  ifelse(t < cp1, int1, ifelse(t < cp2, int1 + time2 * (t - cp1), int3 + time3 * (t - cp2)))
})
pm <- colMeans(mu)

cat("file =", data_file, "\n")
cat("draws_file =", draws_file, "\n")
cat("n =", nrow(dat), "\n")
cat("ndraw =", nrow(draws), "\n")
cat("cp mean checksum =", sprintf("%.15f", mean(cp1) + mean(cp2)), "\n")
cat("param mean checksum =", sprintf("%.15f", sum(colMeans(draws))), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(pm)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(pm[1:10])), "\n")
cat("sigma mean =", sprintf("%.15f", mean(sigma1)), "\n")

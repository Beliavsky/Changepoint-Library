#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmcp_ar_change_data.txt"
draws_file <- "xmcp_ar_change_draws.txt"

source("xmcp_repo_core.R")
env <- load_mcp_repo_core()
mcp <- get("mcp", envir = env)

dat <- read.table(data_file, col.names = c("x", "y"))
model <- list(
  y ~ 1 + ar(2),
  ~ 0 + x + ar(1, 1 + x),
  ~ 0
)

fit <- mcp(
  model = model,
  data = dat,
  sample = "post",
  chains = 3,
  iter = 200,
  adapt = 100,
  cores = 1
)

draws <- do.call(rbind, lapply(fit[["mcmc_post"]], as.matrix))
keep <- c("ar1_1", "ar1_2", "ar1_x_2", "ar2_1", "cp_1", "cp_2", "int_1", "sigma_1", "x_2")
draws <- draws[, keep, drop = FALSE]
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

x <- dat[["x"]]
cp1 <- draws[, "cp_1"]
cp2 <- draws[, "cp_2"]
int1 <- draws[, "int_1"]
x2 <- draws[, "x_2"]
plateau <- int1 + x2 * (cp2 - cp1)

mu <- sapply(x, function(t) {
  ifelse(t < cp1, int1, ifelse(t < cp2, int1 + x2 * (t - cp1), plateau))
})
pm <- colMeans(mu)

cat("file =", data_file, "\n")
cat("draws_file =", draws_file, "\n")
cat("repo =", "C:/rcode/public_domain/github/mcp", "\n")
cat("n =", nrow(dat), "\n")
cat("ndraw =", nrow(draws), "\n")
cat("chains = 3\n")
cat("iter = 200\n")
cat("adapt = 100\n")
cat("cp mean checksum =", sprintf("%.15f", mean(cp1) + mean(cp2)), "\n")
cat("param mean checksum =", sprintf("%.15f", sum(colMeans(draws))), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(pm)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(pm[1:10])), "\n")
cat("ar checksum =", sprintf("%.15f", mean(draws[, "ar1_1"]) + mean(draws[, "ar1_2"]) + mean(draws[, "ar1_x_2"]) + mean(draws[, "ar2_1"])), "\n")
cat("sigma mean =", sprintf("%.15f", mean(draws[, "sigma_1"])), "\n")

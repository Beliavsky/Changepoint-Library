#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmcp_arsigma_data.txt"
draws_file <- "xmcp_arsigma_draws.txt"

source("xmcp_repo_core.R")
env <- load_mcp_repo_core()
mcp <- get("mcp", envir = env)

dat <- read.table(data_file, col.names = c("x", "y"))
model <- list(
  y ~ 1 + ar(1) + sigma(1 + x),
  ~ 0 + x + ar(2, 1 + x) + sigma(1)
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
keep <- c("ar1_1", "ar1_2", "ar1_x_2", "ar2_2", "ar2_x_2", "cp_1", "int_1", "sigma_1", "sigma_2", "sigma_x_1", "x_2")
draws <- draws[, keep, drop = FALSE]
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

x <- dat[["x"]]
cp1 <- draws[, "cp_1"]
int1 <- draws[, "int_1"]
x2 <- draws[, "x_2"]
sigma1 <- draws[, "sigma_1"]
sigma2 <- draws[, "sigma_2"]
sigma_x1 <- draws[, "sigma_x_1"]

mu <- sapply(x, function(t) {
  ifelse(t < cp1, int1, int1 + x2 * (t - cp1))
})
sigma_fit <- sapply(x, function(t) {
  ifelse(t < cp1, pmax(0, sigma1 + sigma_x1 * t), sigma2)
})

pm <- colMeans(mu)
ps <- colMeans(sigma_fit)

cat("file =", data_file, "\n")
cat("draws_file =", draws_file, "\n")
cat("repo =", "C:/rcode/public_domain/github/mcp", "\n")
cat("n =", nrow(dat), "\n")
cat("ndraw =", nrow(draws), "\n")
cat("chains = 3\n")
cat("iter = 200\n")
cat("adapt = 100\n")
cat("cp mean =", sprintf("%.15f", mean(cp1)), "\n")
cat("param mean checksum =", sprintf("%.15f", sum(colMeans(draws))), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(pm)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(pm[1:10])), "\n")
cat("sigma fitted checksum =", sprintf("%.15f", sum(ps)), "\n")
cat("sigma head checksum =", sprintf("%.15f", sum(ps[1:10])), "\n")
cat("ar checksum =", sprintf("%.15f", sum(colMeans(draws[, 1:5, drop = FALSE]))), "\n")
cat("sigma param checksum =", sprintf("%.15f", mean(sigma1) + mean(sigma2) + mean(sigma_x1)), "\n")

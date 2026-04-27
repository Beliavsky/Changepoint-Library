#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmcp_sigma_data.txt"
draws_file <- "xmcp_sigma_draws.txt"

source("xmcp_repo_core.R")
env <- load_mcp_repo_core()
mcp <- get("mcp", envir = env)

model <- list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)

set.seed(30)
dat <- data.frame(x = 1:100)
x <- dat[["x"]]
mu_true <- ifelse(x < 75, 20, 20 + 2 * (x - 75))
sigma_true <- ifelse(x < 25, 7, pmax(0, 25 - 0.45 * (x - 25)))
dat$y <- rnorm(length(x), mu_true, sigma_true)
write.table(dat, file = data_file, row.names = FALSE, col.names = FALSE)

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
keep <- c("cp_1", "cp_2", "int_1", "sigma_1", "sigma_2", "sigma_x_2", "x_3")
draws <- draws[, keep, drop = FALSE]
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

cp1 <- draws[, "cp_1"]
cp2 <- draws[, "cp_2"]
int1 <- draws[, "int_1"]
sigma1 <- draws[, "sigma_1"]
sigma2 <- draws[, "sigma_2"]
sigma_x2 <- draws[, "sigma_x_2"]
x3 <- draws[, "x_3"]

mu <- sapply(x, function(t) {
  ifelse(t < cp1, int1, ifelse(t < cp2, int1, int1 + x3 * (t - cp2)))
})
sigma_fit <- sapply(x, function(t) {
  ifelse(t < cp1, sigma1, pmax(0, sigma2 + sigma_x2 * (t - cp1)))
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
cat("cp mean checksum =", sprintf("%.15f", mean(cp1) + mean(cp2)), "\n")
cat("param mean checksum =", sprintf("%.15f", sum(colMeans(draws))), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(pm)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(pm[1:10])), "\n")
cat("sigma fitted checksum =", sprintf("%.15f", sum(ps)), "\n")
cat("sigma head checksum =", sprintf("%.15f", sum(ps[1:10])), "\n")
cat("sigma param checksum =", sprintf("%.15f", mean(sigma1) + mean(sigma2) + mean(sigma_x2)), "\n")

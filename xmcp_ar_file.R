#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmcp_ar_data.txt"
draws_file <- "xmcp_ar_draws.txt"

source("xmcp_repo_core.R")
env <- load_mcp_repo_core()
mcp <- get("mcp", envir = env)

dat <- read.table(data_file, col.names = c("time", "price"))
fit <- mcp(
  model = list(price ~ 1 + ar(1)),
  data = dat,
  par_x = "time",
  sample = "post",
  chains = 3,
  iter = 200,
  adapt = 100,
  cores = 1
)

draws <- do.call(rbind, lapply(fit[["mcmc_post"]], as.matrix))
keep <- c("ar1_1", "int_1", "sigma_1")
draws <- draws[, keep, drop = FALSE]
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

pm <- rep(mean(draws[, "int_1"]), nrow(dat))

cat("file =", data_file, "\n")
cat("draws_file =", draws_file, "\n")
cat("repo =", "C:/rcode/public_domain/github/mcp", "\n")
cat("n =", nrow(dat), "\n")
cat("ndraw =", nrow(draws), "\n")
cat("chains = 3\n")
cat("iter = 200\n")
cat("adapt = 100\n")
cat("param mean checksum =", sprintf("%.15f", sum(colMeans(draws))), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(pm)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(pm[1:10])), "\n")
cat("ar mean =", sprintf("%.15f", mean(draws[, "ar1_1"])), "\n")
cat("int mean =", sprintf("%.15f", mean(draws[, "int_1"])), "\n")
cat("sigma mean =", sprintf("%.15f", mean(draws[, "sigma_1"])), "\n")

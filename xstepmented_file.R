#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstepmented_data.txt"
draws_file <- "xstepmented_draws.txt"

library(segmented)

dat <- read.table(data_file, col.names = c("x", "y"))
fit0 <- lm(y ~ 1, data = dat)
fit <- stepmented(fit0, seg.Z = ~x, psi = 60)

psi_hat <- as.numeric(fit$psi[, "Est."])
co <- coef(fit)[1:2]
draws <- matrix(c(psi_hat, unname(co)), nrow = 1L)
write.table(draws, file = draws_file, row.names = FALSE, col.names = FALSE)

fv <- fitted(fit)

cat("file =", data_file, "\n")
cat("draws_file =", draws_file, "\n")
cat("n =", nrow(dat), "\n")
cat("psi =", sprintf("%.15f", psi_hat), "\n")
cat("coef checksum =", sprintf("%.15f", sum(co)), "\n")
cat("fitted checksum =", sprintf("%.15f", sum(fv)), "\n")
cat("fitted head checksum =", sprintf("%.15f", sum(fv[1:10])), "\n")
cat("intercept =", sprintf("%.15f", unname(co[1])), "\n")
cat("jump =", sprintf("%.15f", unname(co[2])), "\n")

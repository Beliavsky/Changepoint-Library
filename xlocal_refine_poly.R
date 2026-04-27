#!/usr/bin/env Rscript

data_file <- "xcv_dp_poly_data.txt"
r <- 2L
gamma_set <- c(0.5, 1, 2, 3)
delta <- 5L
delta_lr <- 5L

y <- as.numeric(read.table(data_file, comment.char = "#")[, 1])

t0 <- proc.time()[["elapsed"]]
dp_result <- changepoints::CV.search.DP.poly(y = y, r = r, gamma_set = gamma_set, delta = delta)
min_idx <- which.min(unlist(dp_result$test_error))
cpt_init <- unlist(dp_result$cpt_hat[min_idx])
refined <- changepoints::local.refine.poly(cpt_init, y = y, r = r, delta_lr = delta_lr)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", length(y), "\n")
cat("r =", r, "\n")
cat("gamma_set =", paste(gamma_set, collapse = " "), "\n")
cat("delta =", delta, "\n")
cat("delta_lr =", delta_lr, "\n")
cat("best gamma =", gamma_set[[min_idx]], "\n")
cat("initial changepoints =", paste(cpt_init, collapse = " "), "\n")
cat("refined changepoints =", paste(refined, collapse = " "), "\n")
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

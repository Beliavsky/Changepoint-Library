#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xecp_cp3o_data.txt"
K <- 3L
minsize <- 10L
alpha <- 1

X <- as.matrix(read.table(data_file, comment.char = "#"))
t0 <- proc.time()[["elapsed"]]
fit <- ecp::e.cp3o(X, K = K, minsize = minsize, alpha = alpha, verbose = FALSE)
elapsed <- proc.time()[["elapsed"]] - t0

cat("file =", data_file, "\n")
cat("n =", nrow(X), "\n")
cat("p =", ncol(X), "\n")
cat("K =", K, "\n")
cat("minsize =", minsize, "\n")
cat("alpha =", alpha, "\n")
cat("number =", fit$number, "\n")
cat("estimates =", paste(fit$estimates, collapse = " "), "\n")
cat("gof length =", length(fit$gofM), "\n")
cat(sprintf("gof checksum = %.12f\n", sum(fit$gofM)))
cat(sprintf("cpLoc checksum = %.12f\n", sum(unlist(fit$cpLoc))))
cat("cpLoc lengths =", paste(sapply(fit$cpLoc, length), collapse = " "), "\n")
for (i in seq_along(fit$cpLoc)) {
  cat("cpLoc", i, "=", paste(fit$cpLoc[[i]], collapse = " "), "\n")
}
cat(sprintf("elapsed seconds = %.3f\n", elapsed))

#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmetrics_data.txt"
n_bkps <- 3L
min_size <- 40L
margin <- 40L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
true_bkps <- read_true_bkps(data_file)
n <- length(signal)
width <- max(80L, n %/% 10L)

cat("file =", data_file, "\n")
cat("n =", n, "\n")
cat("true bkps =", paste(true_bkps, collapse = " "), "\n")
cat("n_bkps =", n_bkps, "\n")
cat("min_size =", min_size, "\n")
cat("margin =", margin, "\n")
cat("\n")

methods <- list(
  list(name = "binseg", run = function() binseg_l2_1d(signal, n_bkps, min_size)),
  list(name = "bottomup", run = function() bottomup_l2_1d(signal, n_bkps, min_size, jump = 5L)),
  list(name = "window", run = function() window_l2_1d(signal, n_bkps, width, jump = 5L)),
  list(name = "dynp", run = function() solve_dynp_cost_1d(signal, n_bkps, min_size, "l2"))
)

for (item in methods) {
  t1 <- proc.time()[["elapsed"]]
  bkps <- item$run()
  elapsed <- proc.time()[["elapsed"]] - t1
  pr <- precision_recall_cpt(true_bkps, bkps, margin)
  hd <- hausdorff_cpt(true_bkps, bkps)
  ri <- randindex_cpt(true_bkps, bkps)

  cat("method =", item$name, "\n")
  cat("estimated bkps =", paste(bkps, collapse = " "), "\n")
  cat(sprintf("precision = %.3f\n", pr[["precision"]]))
  cat(sprintf("recall = %.3f\n", pr[["recall"]]))
  cat(sprintf("hausdorff = %.3f\n", hd))
  cat(sprintf("randindex = %.6f\n", ri))
  cat(sprintf("elapsed seconds = %.3f\n", elapsed))
  cat("\n")
}

cat(sprintf("wall time elapsed (s) = %.3f\n", proc.time()[["elapsed"]] - t0))

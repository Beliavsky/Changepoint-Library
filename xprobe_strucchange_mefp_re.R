#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xstrucchange_mefp_re_data.txt"
history_n <- 50L
alpha <- 0.05

dat <- as.data.frame(read.table(data_file, comment.char = "#"))
names(dat) <- c("y", "x")
hist_dat <- dat[1:history_n, , drop = FALSE]

me0 <- strucchange::mefp(y ~ x, type = "RE", data = hist_dat, alpha = alpha)
ans <- strucchange::monitor(me0, data = dat, verbose = FALSE)

proc <- ans[["process"]]
stat <- ans[["statistic"]]
idx <- c(1L, 2L, 10L, 20L, 40L, 70L)

cat(sprintf("critval = %.15f\n", ans[["critval"]]))
cat(sprintf("breakpoint = %d\n", ans[["breakpoint"]]))
for (ii in idx) {
  cat(sprintf("row = %d proc1 = %.15f proc2 = %.15f stat = %.15f\n",
              ii, proc[ii, 1], proc[ii, 2], stat[ii]))
}

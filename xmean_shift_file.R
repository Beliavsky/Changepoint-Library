#!/usr/bin/env Rscript

source("utils.r")

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xmean_shift_data.txt"
pen <- 10.0
jump <- 5L

t0 <- proc.time()[["elapsed"]]
signal <- read_series(data_file)
n <- length(signal)
min_size <- max(5L, n %/% 20L)

prefix <- build_prefix_1d(signal)
cs <- prefix$cs
css <- prefix$css

best_cost <- rep(Inf, n + 1L)
parent <- rep(-1L, n + 1L)
best_cost[1L] <- 0.0
parent[1L] <- 0L

inds <- c(seq.int(0L, n - 1L, by = jump), n)
inds <- inds[inds >= min_size]
admissible <- integer(0)

for (b in inds) {
  new_adm_pt <- ((b - min_size) %/% jump) * jump
  if (new_adm_pt >= 0L && (length(admissible) == 0L || tail(admissible, 1L) != new_adm_pt)) {
    admissible <- c(admissible, new_adm_pt)
  }

  cand_cost <- rep(Inf, length(admissible))
  best_val <- Inf
  best_t <- -1L

  for (idx in seq_along(admissible)) {
    t <- admissible[idx]
    if (!is.finite(best_cost[t + 1L])) next
    seg_cost <- segment_sse_1d(cs, css, t + 1L, b, min_size)
    cand <- best_cost[t + 1L] + seg_cost + pen
    cand_cost[idx] <- cand
    if (cand < best_val) {
      best_val <- cand
      best_t <- t
    }
  }

  best_cost[b + 1L] <- best_val
  parent[b + 1L] <- best_t

  keep <- cand_cost <= best_val + pen
  admissible <- admissible[keep]
}

bkps <- integer(0)
b <- n
while (b > 0L && parent[b + 1L] >= 0L) {
  bkps <- c(b, bkps)
  b <- parent[b + 1L]
}

if (length(bkps) > 0L) {
  bkps <- bkps[-length(bkps)]
}
bkps_out <- c(bkps, n)
elapsed <- proc.time()[["elapsed"]] - t0

print_result_1d(
  data_file,
  n,
  list(
    list(name = "pen", value = pen),
    list(name = "min_size", value = min_size)
  ),
  bkps_out,
  elapsed
)

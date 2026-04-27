#!/usr/bin/env Rscript

perm_cluster_trace <- function(D, points) {
  points <- sort(points)
  n <- nrow(D)
  perm <- seq_len(n)
  K <- length(points) - 1L
  for (i in seq_len(K)) {
    seg <- points[i]:(points[i + 1L] - 1L)
    perm[seg] <- sample(seg)
  }
  list(D = D[perm, perm], perm = perm)
}

sig_test_trace <- function(D, R, changes, min.size, obs) {
  if (R == 0) return(list(p.value = 0, permutations = 0L, perms = NULL))
  over <- 0L
  perms <- matrix(0L, nrow = R, ncol = nrow(D))
  for (f in seq_len(R)) {
    tr <- perm_cluster_trace(D, changes)
    perms[f, ] <- tr$perm
    tmp <- ecp:::e.split(changes, tr$D, min.size, TRUE)
    if (tmp[[4]] >= obs) over <- over + 1L
  }
  list(p.value = (1 + over)/(R + 1), permutations = R, perms = perms)
}

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1L) args[[1L]] else "xecp_edivisive_full_data.txt"
perm_file <- "xecp_edivisive_full_perms.txt"
sig.lvl <- 0.1
R <- 19L
min.size <- 10L
alpha <- 1

X <- as.matrix(read.table(data_file, comment.char = "#"))
n <- nrow(X)
D <- as.matrix(dist(X))^alpha
energy <- new.env(parent = emptyenv())
ecp:::setDim(n, 2, energy)
changes <- c(1L, n + 1L)
k.hat <- 1L
p.values <- permutations <- integer(0)
perm_blocks <- list()
considered.last <- NA_integer_

set.seed(2601)
repeat {
  tmp <- ecp:::e.split(changes, D, min.size, FALSE, energy)
  Estat <- tmp[[4]]
  changes_candidate <- tmp[[3]]
  considered.last <- tail(changes_candidate, 1L)
  if (considered.last == -1L) break
  st <- sig_test_trace(D, R, changes, min.size, Estat)
  p.values <- c(p.values, st$p.value)
  permutations <- c(permutations, st$permutations)
  perm_blocks[[length(perm_blocks) + 1L]] <- st$perms
  if (st$p.value > sig.lvl) break
  changes <- changes_candidate
  k.hat <- k.hat + 1L
}

perm_mat <- if (length(perm_blocks) > 0L) do.call(rbind, perm_blocks) else matrix(integer(0), nrow = 0L, ncol = n)
write.table(perm_mat, file = perm_file, row.names = FALSE, col.names = FALSE)

estimates <- sort(changes)
cluster <- rep(seq_along(diff(estimates)), diff(estimates))

cat("file =", data_file, "\n")
cat("perm_file =", perm_file, "\n")
cat("n =", nrow(X), "\n")
cat("p =", ncol(X), "\n")
cat("sig.lvl =", sig.lvl, "\n")
cat("R =", R, "\n")
cat("min.size =", min.size, "\n")
cat("alpha =", alpha, "\n")
cat("k.hat =", k.hat, "\n")
cat("estimates =", paste(estimates, collapse = " "), "\n")
cat("order.found =", paste(changes, collapse = " "), "\n")
cat("considered.last =", considered.last, "\n")
cat("p.values =", paste(sprintf("%.12f", p.values), collapse = " "), "\n")
cat("permutations =", paste(permutations, collapse = " "), "\n")
cat(sprintf("cluster checksum = %.12f\n", sum(cluster)))

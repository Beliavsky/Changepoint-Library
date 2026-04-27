#' Read a numeric univariate series from a text file, skipping comment lines.
read_series <- function(path) {
  vals <- scan(path, what = numeric(), comment.char = "#", quiet = TRUE)
  if (length(vals) == 0L) stop("input file contains no numeric data")
  vals
}

#' Read true breakpoint metadata from a comment line '# true_bkps = ...'.
read_true_bkps <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hit <- grep("^#\\s*true_bkps\\s*=", lines, value = TRUE)
  if (length(hit) == 0L) stop("input file is missing '# true_bkps =' metadata")
  as.integer(strsplit(sub("^#\\s*true_bkps\\s*=\\s*", "", hit[1L]), "\\s+")[[1L]])
}

#' Read a numeric matrix from a text file with comment lines starting with '#'.
read_matrix <- function(path) {
  x <- as.matrix(read.table(path, comment.char = "#"))
  if (length(x) == 0L) stop("input file contains no numeric data")
  x
}

#' Build 1D cumulative sums and cumulative squared sums for fast segment costs.
build_prefix_1d <- function(x) {
  list(
    cs = c(0, cumsum(x)),
    css = c(0, cumsum(x * x))
  )
}

#' Build multivariate cumulative sums and cumulative squared sums by column.
build_prefix_mv <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  sx <- rbind(rep(0, p), apply(x, 2, cumsum))
  sxx <- rbind(rep(0, p), apply(x * x, 2, cumsum))
  list(sx = sx, sxx = sxx)
}

#' Compute univariate l2 segment SSE from prefix sums over inclusive indices.
segment_sse_1d <- function(cs, css, start_idx, end_idx, min_size) {
  len <- end_idx - start_idx + 1L
  if (len < min_size) return(Inf)
  sum_z <- cs[end_idx + 1L] - cs[start_idx]
  sum_zz <- css[end_idx + 1L] - css[start_idx]
  sum_zz - (sum_z * sum_z) / len
}

#' Compute multivariate l2 segment SSE from prefix sums over inclusive indices.
segment_sse_mv <- function(sx, sxx, start_idx, end_idx, min_size) {
  len <- end_idx - start_idx + 1L
  if (len < min_size) return(Inf)
  sum_x <- sx[end_idx + 1L, , drop = FALSE] - sx[start_idx, , drop = FALSE]
  sum_xx <- sxx[end_idx + 1L, , drop = FALSE] - sxx[start_idx, , drop = FALSE]
  sum(sum_xx - (sum_x * sum_x) / len)
}

#' Compute univariate l2 segment SSE on a half-open interval [start0, end0).
segment_sse_half_open_1d <- function(cs, css, start0, end0) {
  len <- end0 - start0
  if (len <= 0L) return(0.0)
  sum_z <- cs[end0 + 1L] - cs[start0 + 1L]
  sum_zz <- css[end0 + 1L] - css[start0 + 1L]
  sum_zz - (sum_z * sum_z) / len
}

#' Compute the univariate normal cost on inclusive indices.
segment_normal_cost_1d <- function(cs, css, start_idx, end_idx, min_size, small_diag = 1e-6) {
  len <- end_idx - start_idx + 1L
  if (len < min_size) return(Inf)
  sum_z <- cs[end_idx + 1L] - cs[start_idx]
  sum_zz <- css[end_idx + 1L] - css[start_idx]
  var_hat <- max(sum_zz / len - (sum_z / len)^2, 0.0)
  len * log(var_hat + small_diag)
}

#' Compute the univariate l1 cost on inclusive indices.
segment_l1_cost_1d <- function(signal, start_idx, end_idx, min_size) {
  len <- end_idx - start_idx + 1L
  if (len < min_size) return(Inf)
  seg <- signal[start_idx:end_idx]
  med <- median(seg)
  sum(abs(seg - med))
}

#' Solve exact fixed-K dynamic programming for one univariate cost model.
solve_dynp_cost_1d <- function(signal, n_bkps, min_size, model) {
  n <- length(signal)
  max_m <- n_bkps + 1L
  prefix <- build_prefix_1d(signal)
  cs <- prefix$cs
  css <- prefix$css

  if (model == "rbf") {
    pref <- build_rbf_prefix_mv(matrix(signal, ncol = 1L), rbf_gamma_from_matrix(matrix(signal, ncol = 1L)))
  } else {
    pref <- NULL
  }

  seg_cost <- function(start_idx, end_idx) {
    if (model == "l1") {
      segment_l1_cost_1d(signal, start_idx, end_idx, min_size)
    } else if (model == "l2") {
      segment_sse_1d(cs, css, start_idx, end_idx, min_size)
    } else if (model == "normal") {
      segment_normal_cost_1d(cs, css, start_idx, end_idx, min_size)
    } else if (model == "rbf") {
      segment_rbf_cost_mv(pref, start_idx - 1L, end_idx, min_size)
    } else {
      stop("unsupported model")
    }
  }

  dp <- matrix(Inf, nrow = n, ncol = max_m)
  parent <- matrix(0L, nrow = n, ncol = max_m)

  for (i in seq_len(n)) {
    dp[i, 1L] <- seg_cost(1L, i)
  }

  if (max_m >= 2L) {
    for (m in 2:max_m) {
      for (i in m:n) {
        if (i < m * min_size) next
        best_val <- Inf
        best_k <- 0L
        for (k in (m - 1L):(i - 1L)) {
          if (!is.finite(dp[k, m - 1L])) next
          cand <- dp[k, m - 1L] + seg_cost(k + 1L, i)
          if (cand < best_val) {
            best_val <- cand
            best_k <- k
          }
        }
        dp[i, m] <- best_val
        parent[i, m] <- best_k
      }
    }
  }

  seg_ends <- integer(max_m)
  idx <- n
  seg_ends[max_m] <- n
  if (max_m >= 2L) {
    for (m in seq(max_m, 2L, by = -1L)) {
      seg_ends[m - 1L] <- parent[idx, m]
      idx <- seg_ends[m - 1L]
    }
  }
  seg_ends
}

#' Run binary segmentation with l2 cost on a univariate signal.
binseg_l2_1d <- function(signal, n_bkps, min_size) {
  prefix <- build_prefix_1d(signal)
  cs <- prefix$cs
  css <- prefix$css

  seg_gain <- function(start0, end0) {
    best_gain <- -Inf
    best_k <- NA_integer_
    full_cost <- segment_sse_half_open_1d(cs, css, start0, end0)
    lo <- start0 + min_size
    hi <- end0 - min_size
    if (lo > hi) return(c(bkp = NA_integer_, gain = -Inf))
    for (k in seq.int(lo, hi)) {
      gain <- full_cost - segment_sse_half_open_1d(cs, css, start0, k) -
        segment_sse_half_open_1d(cs, css, k, end0)
      if (gain > best_gain) {
        best_gain <- gain
        best_k <- k
      }
    }
    c(bkp = best_k, gain = best_gain)
  }

  segs <- list(c(0L, length(signal)))
  bkps <- integer(0)
  for (step in seq_len(n_bkps)) {
    best <- list(seg = NA_integer_, bkp = NA_integer_, gain = -Inf)
    for (i in seq_along(segs)) {
      cand <- seg_gain(segs[[i]][1L], segs[[i]][2L])
      if (is.na(cand[["bkp"]])) next
      if (cand[["gain"]] > best$gain) {
        best <- list(seg = i, bkp = as.integer(cand[["bkp"]]), gain = cand[["gain"]])
      }
    }
    if (is.na(best$bkp)) break
    seg <- segs[[best$seg]]
    segs <- append(segs[-best$seg], list(c(seg[1L], best$bkp), c(best$bkp, seg[2L])))
    bkps <- sort(c(bkps, best$bkp))
  }
  c(bkps, length(signal))
}

#' Refine internal breakpoints by exact local l2 re-optimization.
refine_bkps_l2_1d <- function(signal, bkps, min_size, max_iter = 20L) {
  if (length(bkps) == 0L) return(integer(0))

  n <- length(signal)
  if (tail(bkps, 1L) != n) stop("bkps must include the final endpoint n")

  prefix <- build_prefix_1d(signal)
  cs <- prefix$cs
  css <- prefix$css
  refined <- as.integer(bkps)

  for (iter in seq_len(max_iter)) {
    changed <- FALSE
    for (i in seq_len(length(refined) - 1L)) {
      left0 <- if (i == 1L) 0L else refined[i - 1L]
      right0 <- refined[i + 1L]
      lo <- left0 + min_size
      hi <- right0 - min_size
      if (lo > hi) next

      best_k <- refined[i]
      best_cost <- segment_sse_half_open_1d(cs, css, left0, best_k) +
        segment_sse_half_open_1d(cs, css, best_k, right0)

      for (k in seq.int(lo, hi)) {
        cand_cost <- segment_sse_half_open_1d(cs, css, left0, k) +
          segment_sse_half_open_1d(cs, css, k, right0)
        if (cand_cost < best_cost) {
          best_cost <- cand_cost
          best_k <- k
        }
      }

      if (best_k != refined[i]) {
        refined[i] <- best_k
        changed <- TRUE
      }
    }
    if (!changed) break
  }

  refined
}

#' Build the initial regular binary partition used by ruptures BottomUp.
grow_bottomup_leaves_1d <- function(start0, end0, min_size, jump) {
  mid <- 0.5 * (start0 + end0)
  bkps <- integer(0)
  for (bkp in seq.int(start0, end0 - 1L)) {
    if (bkp %% jump != 0L) next
    if (bkp - start0 < min_size) next
    if (end0 - bkp < min_size) next
    bkps <- c(bkps, bkp)
  }

  if (length(bkps) == 0L) {
    return(list(c(start0, end0)))
  }

  bkp <- bkps[which.min(abs(bkps - mid))]
  c(
    grow_bottomup_leaves_1d(start0, bkp, min_size, jump),
    grow_bottomup_leaves_1d(bkp, end0, min_size, jump)
  )
}

#' Run bottom-up l2 segmentation using the same grow-then-merge structure as ruptures.
bottomup_l2_1d <- function(signal, n_bkps, min_size, jump = 5L) {
  n <- length(signal)
  prefix <- build_prefix_1d(signal)
  cs <- prefix$cs
  css <- prefix$css

  leaves <- grow_bottomup_leaves_1d(0L, n, min_size, jump)
  seg_start <- vapply(leaves, `[[`, integer(1), 1L)
  seg_end <- vapply(leaves, `[[`, integer(1), 2L)
  seg_cost <- mapply(
    function(a, b) segment_sse_half_open_1d(cs, css, a, b),
    seg_start,
    seg_end
  )

  nseg <- length(seg_end)
  while (nseg > n_bkps + 1L) {
    best_gain <- Inf
    best_i <- NA_integer_
    best_cost <- Inf

    for (i in seq_len(nseg - 1L)) {
      merged_cost <- segment_sse_half_open_1d(cs, css, seg_start[i], seg_end[i + 1L])
      gain <- merged_cost - seg_cost[i] - seg_cost[i + 1L]
      if (gain < best_gain) {
        best_gain <- gain
        best_i <- i
        best_cost <- merged_cost
      }
    }

    seg_end[best_i] <- seg_end[best_i + 1L]
    seg_cost[best_i] <- best_cost
    if (best_i < nseg - 1L) {
      keep <- seq_len(nseg) != (best_i + 1L)
      seg_start <- seg_start[keep]
      seg_end <- seg_end[keep]
      seg_cost <- seg_cost[keep]
    } else {
      seg_start <- seg_start[-(best_i + 1L)]
      seg_end <- seg_end[-(best_i + 1L)]
      seg_cost <- seg_cost[-(best_i + 1L)]
    }
    nseg <- nseg - 1L
  }

  as.integer(seg_end)
}

#' Run the sliding-window l2 method on a univariate signal.
window_l2_1d <- function(signal, n_bkps, width, jump = 5L) {
  n <- length(signal)
  width <- 2L * (width %/% 2L)
  prefix <- build_prefix_1d(signal)
  cs <- prefix$cs
  css <- prefix$css

  inds <- seq.int(0L, n - 1L, by = jump)
  keep <- inds >= width %/% 2L & inds < n - width %/% 2L
  inds <- inds[keep]

  if (length(inds) == 0L) return(as.integer(n))

  score <- numeric(length(inds))
  for (i in seq_along(inds)) {
    k <- inds[i]
    start0 <- k - width %/% 2L
    end0 <- k + width %/% 2L
    full_cost <- segment_sse_1d(cs, css, start0 + 1L, end0, 1L)
    left_cost <- segment_sse_1d(cs, css, start0 + 1L, k, 1L)
    right_cost <- segment_sse_1d(cs, css, k + 1L, end0, 1L)
    score[i] <- full_cost - left_cost - right_cost
  }

  order <- max(width %/% (2L * jump), 1L)
  peaks <- find_window_peaks(score, inds, order)
  if (length(peaks$bkps) == 0L) return(as.integer(n))
  ord <- order(peaks$gain)
  peak_bkps <- peaks$bkps[ord]
  chosen <- tail(peak_bkps, min(n_bkps, length(peak_bkps)))
  c(sort(chosen), n)
}

#' Compute precision and recall for breakpoint lists with a tolerance margin.
precision_recall_cpt <- function(true_bkps, est_bkps, margin) {
  if (length(est_bkps) == 1L) return(c(precision = 0, recall = 0))
  used <- rep(FALSE, length(est_bkps) - 1L)
  tp <- 0L
  for (tb in true_bkps[-length(true_bkps)]) {
    for (i in seq_len(length(est_bkps) - 1L)) {
      eb <- est_bkps[i]
      if (!used[i] && eb - margin < tb && tb < eb + margin) {
        used[i] <- TRUE
        tp <- tp + 1L
        break
      }
    }
  }
  c(
    precision = tp / (length(est_bkps) - 1L),
    recall = tp / (length(true_bkps) - 1L)
  )
}

#' Compute the Hausdorff distance between two breakpoint sets.
hausdorff_cpt <- function(bkps1, bkps2) {
  a <- bkps1[-length(bkps1)]
  b <- bkps2[-length(bkps2)]
  if (length(a) == 0L || length(b) == 0L) return(Inf)
  d <- abs(outer(a, b, "-"))
  max(max(apply(d, 1L, min)), max(apply(d, 2L, min)))
}

#' Compute the Rand index between two segmentations.
randindex_cpt <- function(bkps1, bkps2) {
  n_samples <- tail(bkps1, 1L)
  bkps1_with_0 <- c(0L, bkps1)
  bkps2_with_0 <- c(0L, bkps2)
  disagreement <- 0
  beginj <- 1L

  for (i in seq_len(length(bkps1))) {
    start1 <- bkps1_with_0[i]
    end1 <- bkps1_with_0[i + 1L]
    for (j in beginj:length(bkps2)) {
      start2 <- bkps2_with_0[j]
      end2 <- bkps2_with_0[j + 1L]
      nij <- max(min(end1, end2) - max(start1, start2), 0L)
      disagreement <- disagreement + nij * abs(end1 - end2)
      if (end1 < end2) {
        break
      } else {
        beginj <- j + 1L
      }
    }
  }

  1.0 - disagreement / (n_samples * (n_samples - 1) / 2)
}

#' Print a standard result block for 1D comparison scripts.
print_result_1d <- function(data_file, n, extras, bkps, elapsed) {
  cat("file =", data_file, "\n")
  cat("n =", n, "\n")
  for (item in extras) {
    cat(item$name, "=", item$value, "\n")
  }
  cat("estimated bkps =", paste(bkps, collapse = " "), "\n")
  cat(sprintf("elapsed seconds = %.3f\n", elapsed))
}

#' Print a standard result block for multivariate comparison scripts.
print_result_mv <- function(data_file, n, p, extras, bkps, elapsed) {
  cat("file =", data_file, "\n")
  cat("n =", n, "\n")
  cat("p =", p, "\n")
  for (item in extras) {
    cat(item$name, "=", item$value, "\n")
  }
  cat("estimated bkps =", paste(bkps, collapse = " "), "\n")
  cat(sprintf("elapsed seconds = %.3f\n", elapsed))
}

#' Estimate the RBF gamma parameter from the median pairwise squared distance.
rbf_gamma_from_matrix <- function(x) {
  d <- dist(x, method = "euclidean")
  d2 <- as.numeric(d)^2
  med <- median(d2)
  if (!is.finite(med) || med == 0) return(1.0)
  1.0 / med
}

#' Build a 2D prefix sum of the RBF Gram matrix for a multivariate signal.
build_rbf_prefix_mv <- function(x, gamma) {
  d2 <- as.matrix(dist(x, method = "euclidean"))^2
  scaled <- pmin(pmax(gamma * d2, 1e-2), 1e2)
  gram <- exp(-scaled)
  diag(gram) <- 1.0
  pref <- matrix(0, nrow = nrow(gram) + 1L, ncol = ncol(gram) + 1L)
  pref[-1L, -1L] <- apply(apply(gram, 2L, cumsum), 1L, cumsum)
  t(pref)
}

#' Compute the RBF-kernel segment cost on a half-open interval [start0, end0).
segment_rbf_cost_mv <- function(pref, start0, end0, min_size) {
  len <- end0 - start0
  if (len < min_size) return(Inf)
  s <- start0 + 1L
  e <- end0
  sub_sum <- pref[e + 1L, e + 1L] - pref[s, e + 1L] - pref[e + 1L, s] + pref[s, s]
  len - sub_sum / len
}

#' Find strict local maxima with wraparound using a symmetric window order.
find_window_peaks <- function(score, inds, order) {
  n <- length(score)
  if (n == 0L) {
    return(list(bkps = integer(0), gain = numeric(0)))
  }

  is_peak <- logical(n)
  for (i in seq_len(n)) {
    peak <- TRUE
    for (j in seq.int(-order, order)) {
      if (j == 0L) next
      idx <- ((i - 1L + j) %% n) + 1L
      if (score[i] <= score[idx]) {
        peak <- FALSE
        break
      }
    }
    is_peak[i] <- peak
  }

  list(bkps = inds[is_peak], gain = score[is_peak])
}

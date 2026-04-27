library(changepoints)
library(gglasso)

probe_case <- function(DATA, s = 20L, e = 60L, eta_set = c(25L, 38L, 45L, 48L), zeta_set = c(0.01, 0.05)) {
  DATA.temp <- DATA
  if (ncol(DATA) %% 2 == 0) {
    DATA.temp <- DATA[, 2:ncol(DATA)]
  }
  N <- ncol(DATA.temp)
  p <- nrow(DATA.temp)
  X_curr <- DATA.temp[, 1:(N - 1), drop = FALSE]
  X_futu <- DATA.temp[, 2:N, drop = FALSE]
  X_curr.train <- X_curr[, seq(1, N - 1, 2), drop = FALSE]
  X_curr.test <- X_curr[, seq(2, N - 1, 2), drop = FALSE]
  X_futu.train <- X_futu[, seq(1, N - 1, 2), drop = FALSE]
  X_futu.test <- X_futu[, seq(2, N - 1, 2), drop = FALSE]

  cat("probe s            =", s, "\n")
  cat("probe e            =", e, "\n")
  for (eta in eta_set) {
    eta.rel <- eta - s + 2L
    for (zeta in zeta_set) {
      cat("probe mode         = objective\n")
      cat("probe eta          =", eta, "\n")
      cat("probe zeta         =", format(zeta, digits = 10), "\n")
      X.convert <- changepoints:::X.glasso.converter.VAR1(X_curr.train[, s:e, drop = FALSE], eta, s)
      for (m in seq_len(p)) {
        y.vec <- X_futu.train[m, s:e]
        fit <- gglasso(x = X.convert, y = y.vec, group = rep(1:p, 2), intercept = FALSE, loss = "ls", lambda = zeta / ncol(X_curr.train[, s:e, drop = FALSE]))
        beta <- as.vector(fit$beta)
        obj.term <- sum((y.vec - X.convert %*% beta)^2) + zeta * sqrt(sum(beta^2))
        cat("probe response     =", m, "\n")
        cat("probe beta         =", paste(format(beta, digits = 16), collapse = " "), "\n")
        cat("probe obj term     =", format(obj.term, digits = 16), "\n")
      }

      cat("probe mode         = test\n")
      cat("probe eta          =", eta, "\n")
      cat("probe eta rel      =", eta.rel, "\n")
      cat("probe zeta         =", format(zeta, digits = 10), "\n")
      X.convert <- changepoints:::X.glasso.converter.VAR1(X_curr.train[, s:e, drop = FALSE], eta.rel, 1)
      X.test.convert <- changepoints:::X.glasso.converter.VAR1(X_curr.test[, s:e, drop = FALSE], eta.rel, 1)
      for (m in seq_len(p)) {
        y.train <- X_futu.train[m, s:e]
        y.test <- X_futu.test[m, s:e]
        fit <- gglasso(x = X.convert, y = y.train, group = rep(1:p, each = 2), intercept = FALSE, loss = "ls", lambda = zeta / ncol(X_curr.train[, s:e, drop = FALSE]))
        beta <- as.vector(fit$beta)
        test.term <- sum((y.test - X.test.convert %*% beta)^2) + zeta * sqrt(sum(beta^2))
        cat("probe response     =", m, "\n")
        cat("probe beta         =", paste(format(beta, digits = 16), collapse = " "), "\n")
        cat("probe test term    =", format(test.term, digits = 16), "\n")
      }
    }
  }
}

source_data <- function(path = "xcv_dp_var1_data.txt") {
  as.matrix(read.table(path, comment.char = "#"))
}

args <- commandArgs(trailingOnly = TRUE)
data_file <- if (length(args) >= 1) args[1] else "xcv_dp_var1_data.txt"
DATA <- t(source_data(data_file))
cat("file               =", data_file, "\n")
cat("p                  =", nrow(DATA), "\n")
cat("n                  =", ncol(DATA), "\n")
probe_case(DATA)

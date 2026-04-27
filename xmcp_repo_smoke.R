args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "C:/rcode/public_domain/github/mcp"

source("xmcp_repo_core.R")
env <- load_mcp_repo_core(repo)
mcp <- get("mcp", envir = env)

set.seed(1)
data <- data.frame(
  x = 1:20,
  y = c(rep(2, 7), seq(2, 8, length.out = 6), seq(10, 14, length.out = 7)) +
    rnorm(20, sd = 0.2)
)
model <- list(
  y ~ 1,
  ~ 0 + x,
  ~ 1 + x
)

fit <- mcp(
  model,
  data = data,
  sample = "none"
)

cat("repo =", normalizePath(repo, winslash = "/"), "\n")
cat("jags status =", mcp_repo_jags_status(), "\n")
cat("fit class =", paste(class(fit), collapse = " "), "\n")
cat("has jags_code =", !is.null(fit$jags_code), "\n")
cat("jags_code checksum =", sum(utf8ToInt(paste(fit$jags_code, collapse = "\n"))), "\n")

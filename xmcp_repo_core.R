load_mcp_repo_core <- function(repo = "C:/rcode/public_domain/github/mcp") {
  rdir <- file.path(repo, "R")
  files <- c(
    "misc.R",
    "is_assert.R",
    "families.R",
    "data.R",
    "comparison.R",
    "get_formula.R",
    "get_segment_table.R",
    "get_prior.R",
    "get_jagscode.R",
    "run_jags.R",
    "mcpfit_methods.R",
    "mcp.R"
  )

  missing <- files[!file.exists(file.path(rdir, files))]
  if (length(missing) > 0) {
    stop("Missing mcp repo files: ", paste(missing, collapse = ", "))
  }

  libs <- c(
    "magrittr",
    "dplyr",
    "rlang",
    "stringr",
    "tibble",
    "tidyr",
    "tidyselect",
    "future",
    "future.apply",
    "coda"
  )
  for (pkg in libs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  env <- new.env(parent = globalenv())
  for (f in files) {
    sys.source(file.path(rdir, f), envir = env)
  }

  env
}

mcp_repo_jags_status <- function() {
  tryCatch({
    suppressPackageStartupMessages(library(rjags))
    "ok"
  }, error = function(e) {
    conditionMessage(e)
  })
}


suppressPackageStartupMessages({
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  if (!require("kmdata", quietly = TRUE)) {
    cat("Installing kmdata from GitHub (raredd/kmdata)...\n")
    remotes::install_github("raredd/kmdata", quiet = TRUE)
  }
  library(kmdata)
})

trials  <- commandArgs(trailingOnly = TRUE)
out_dir <- trials[length(trials)]   # last arg is output dir
trials  <- trials[-length(trials)]  # remaining args are trial names

all_names <- ls("package:kmdata")
cat(sprintf("kmdata contains %d datasets.\n", length(all_names)))
cat(sprintf("Requested: %s\n", paste(trials, collapse = ", ")))

for (trial in trials) {
  # Exact or fuzzy match against package objects
  matched <- trial
  if (!matched %in% all_names) {
    pat <- gsub("_", ".", trial, fixed = TRUE)
    fuzzy <- all_names[grepl(pat, all_names, ignore.case = TRUE)]
    if (length(fuzzy) == 0) {
      cat(sprintf("[SKIP] '%s' not found. Fuzzy search returned nothing.\n", trial))
      next
    }
    matched <- fuzzy[1]
    cat(sprintf("[INFO] Matched '%s' -> '%s'.\n", trial, matched))
  }
  df <- get(matched, envir = asNamespace("kmdata"))
  if (!is.data.frame(df)) {
    cat(sprintf("[SKIP] '%s' is not a data frame.\n", matched))
    next
  }
  fn <- file.path(out_dir, paste0(trial, ".csv"))
  write.csv(df, fn, row.names = FALSE)
  cat(sprintf("[OK]   %s -> %s  (%d rows, cols: %s)\n",
              matched, fn, nrow(df), paste(colnames(df), collapse = ", ")))
}

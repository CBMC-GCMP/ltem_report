# Utility functions for LTEM pipeline organization, flags, saving, and logging

safe_mkdir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Create required directory structure, idempotently
make_dirs <- function(cfg) {
  out <- cfg$outputs_dir
  dirs <- c(
    file.path(out, "temp"),
    file.path(out, "temp", "species_names"),
    file.path(out, "temp", "reefs"),
    file.path(out, "review"),
    file.path(out, "review", "size_check"),
    file.path(out, "review", "meta_check"),
    file.path(out, "review", "transect_check")
  )
  invisible(lapply(dirs, safe_mkdir))
}

.timestamp_str <- function() format(Sys.time(), "%Y%m%d-%H%M")

.stage_filename <- function(nn, slug, fmt) {
  sprintf("stage-%02d_%s_%s.%s", as.integer(nn), slug, .timestamp_str(), tolower(fmt))
}

# Save a stage artifact in a standard way
# dir_path: directory to write into
# nn: integer stage number
# slug: short slug for the kind of artifact
# fmt: one of 'rds', 'csv', 'parquet'
# returns the full path written
save_stage <- function(df, dir_path, nn, slug, fmt = "rds") {
  safe_mkdir(dir_path)
  fname <- .stage_filename(nn, slug, fmt)
  fpath <- file.path(dir_path, fname)
  fmt <- tolower(fmt)
  if (fmt == "rds") {
    saveRDS(df, fpath)
  } else if (fmt == "csv") {
    utils::write.csv(df, fpath, row.names = FALSE)
  } else if (fmt == "parquet") {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      warning("arrow not installed; writing RDS instead of parquet: ", fpath)
      fpath <- sub("\\.parquet$", ".rds", fpath)
      saveRDS(df, fpath)
    } else {
      arrow::write_parquet(df, fpath)
    }
  } else {
    stop("Unsupported format: ", fmt)
  }
  return(fpath)
}

# Append a flag token to rows indicated by logical idx
append_flag <- function(df, idx, token) {
  if (is.null(df$Flag)) df$Flag <- ""
  token <- as.character(token)
  if (!is.logical(idx) || length(idx) != nrow(df)) {
    stop("idx must be a logical vector of length nrow(df)")
  }
  add_token <- function(x) {
    x <- as.character(x %||% "")
    if (nzchar(x)) paste0(x, ";", token) else token
  }
  df$Flag[idx] <- vapply(df$Flag[idx], add_token, character(1))
  df
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Summarize flags: returns named integer vector with counts per token
summarize_flags <- function(df) {
  if (is.null(df$Flag) || all(!nzchar(df$Flag))) return(integer(0))
  toks <- unlist(strsplit(df$Flag, ";", fixed = TRUE))
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0) return(integer(0))
  sort(table(toks), decreasing = TRUE)
}

# Bilingual stage completion logger
stage_log_en_es <- function(nn, en_title, es_title, out_dir, rows_updated, flag_counts, next_hint_en, next_hint_es) {
  flags_str <- if (length(flag_counts)) paste(sprintf("%s (%d)", names(flag_counts), as.integer(flag_counts)), collapse = ", ") else "None"
  msg <- paste0(
    "\n",
    "✅ EN: Stage ", nn, " complete — ", en_title, ".\n",
    "📁 Outputs: ", out_dir, "\n",
    "🧾 Changes: ", rows_updated, " rows updated\n",
    "🏷️ Flags: ", flags_str, "\n",
    "➡️ Next: ", next_hint_en, "\n\n",
    "✅ ES: Etapa ", nn, " completada — ", es_title, ".\n",
    "📁 Salidas: ", out_dir, "\n",
    "🧾 Cambios: ", rows_updated, " filas actualizadas\n",
    "🏷️ Banderas: ", flags_str, "\n",
    "➡️ Siguiente: ", next_hint_es, "\n"
  )
  message(msg)
}

# Write a sidecar .txt explaining a review output
write_sidecar <- function(out_file, text) {
  info_path <- paste0(tools::file_path_sans_ext(out_file), ".txt")
  cat(text, file = info_path)
  return(info_path)
}

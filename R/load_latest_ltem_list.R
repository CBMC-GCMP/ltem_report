library(tidyverse)
library(readxl)
library(readr)
library(lubridate)

## Load the most recent LTEM list file
#'
#' Selects and loads the most recent LTEM list file (e.g., species, reefs,
#' sizes, DOV reefs) from a directory. "Most recent" is determined first by
#' parsing a date from the filename (e.g., "2025-05-06" or "2025_04-21"); if no
#' valid date is found, file modification time is used instead.
#'
#' The function supports both `.xlsx` and `.csv` files and chooses the reader
#' accordingly.
#'
#' @param keyword Character scalar used to match list files by name
#'   (e.g., "species", "reef", "size", "dovs"). Matching is case-insensitive.
#' @param path Character scalar, directory where list files are stored.
#'   Default is `"data/lists/"`.
#' @param date_regex Regular expression used to extract a date from the
#'   filename. Default matches patterns like `"2025-05-06"` or `"2025_05-06"`.
#'
#' @return A data.frame/tibble with the contents of the most recent matching
#'   list file.
#'
#' @details
#' - Files are first filtered by `keyword` (case-insensitive) within `path`.
#' - If the filename contains a date matching `date_regex`, this date is parsed
#'   (after normalizing `_`/`-` separators) and used to identify the most recent
#'   file.
#' - If no valid date is found in any matching filenames, the function falls
#'   back to file modification time (`file.info()$mtime`).
#' - `.xlsx` files are loaded with [readxl::read_xlsx()], `.csv` files with
#'   [readr::read_csv()].
#'
#' @examples
#' # Load most recent species catalog
#' # species_list <- load_latest_ltem_list("species")
#'
#' # Load most recent reefs catalog
#' # reefs_list <- load_latest_ltem_list("reef")
#'
#' # Load most recent size list
#' # size_list <- load_latest_ltem_list("size")
#'
#' # Load most recent DOV reefs list
#' # dovs_list <- load_latest_ltem_list("dovs")
#'
#' @keywords data-import files LTEM
load_latest_ltem_list <- function(
    keyword,
    path = "data/lists/updates",
    date_regex = "\\d{4}[-_]\\d{2}[-_]\\d{2}"  # matches 2025-05-06 or 2025_04-21
) {
  # List files containing the keyword (case-insensitive)
  files <- list.files(
    path,
    pattern = keyword,
    full.names = TRUE,
    ignore.case = TRUE
  )
  files <- files[!startsWith(basename(files), "~$")]
  
  if (length(files) == 0) {
    stop("No files found matching keyword '", keyword, "' in path ", path)
  }
  
  # Build a table with file metadata
  file_tbl <- tibble(file = files) %>%
    mutate(
      base_name      = basename(file),
      date_str       = str_extract(base_name, date_regex),
      date_str_clean = str_replace_all(date_str, "_", "-"),
      date           = suppressWarnings(ymd(date_str_clean)),
      modified       = file.info(file)$mtime,
      ext            = tools::file_ext(file)
    )
  
  # Choose most recent file by date if available, otherwise by modification time
  if (any(!is.na(file_tbl$date))) {
    latest_file <- file_tbl %>%
      filter(!is.na(date)) %>%
      arrange(desc(date)) %>%
      slice(1) %>%
      pull(file)
  } else {
    latest_file <- file_tbl %>%
      arrange(desc(modified)) %>%
      slice(1) %>%
      pull(file)
  }
  
  message("Loading LTEM list from: ", latest_file)
  
  # Dispatch to appropriate reader based on extension
  if (grepl("\\.xlsx$", latest_file, ignore.case = TRUE)) {
    out <- read_xlsx(latest_file)
  } else if (grepl("\\.csv$", latest_file, ignore.case = TRUE)) {
    out <- read_csv(latest_file, show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", latest_file)
  }
  
  return(out)
}

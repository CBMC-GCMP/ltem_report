suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyselect)
})
#' Build a flag table of modified records
#'
#' Filters records with status containing "Modified_" and selects identifiers
#' and correction columns to produce a concise report. Prints unique values
#' of `correct_sp` for quick inspection.
#'
#' @param LTEM_corrected Data frame with corrected LTEM data, including
#'   columns `Status`, identifiers (e.g., `IDReef`, `IDSize`, `IDSpecies`),
#'   context columns (e.g., `Label`, `Region`, `Transect`, etc.), and any
#'   `*_Before`/`correct_*` columns.
#' @return A data frame of flagged rows and selected columns.
#' @examples
#' # flagged <- flags(ltem_corrected)
flags <- function(LTEM_corrected) {
  IDs <- c("IDReef", "IDSize", "IDSpecies")
  correct <- c("correct_size", "correct_reef", "correct_id", "correct_sp")
  Before <- c("IDBefore", "ReefBefore", "SpeciesBefore")
  flags <- LTEM_corrected %>%
    filter(str_detect(Status, "Modified_")) %>%
    select(
      Label,
      Year,
      Month,
      Day,
      Region,
      Depth,
      Transect,
      Observer,
      Habitat,
      any_of(IDs),
      Species,
      any_of(Before),
      any_of(correct),
      Status,
      # SpeciesBefore, correct_sp
    )
  print(unique(flags$correct_sp))
  
  output <- flags
  
  return(output)
}

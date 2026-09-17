library(tidyverse)
library(writexl)
library(lubridate)

## Save an updated LTEM list with an auto-generated date-stamped filename
#'
#' Saves an LTEM list (species, reefs, sizes, dovs reefs, etc.) to the updates
#' directory with a standardized filename containing the current date.
#'
#' @param data A data.frame/tibble to save.
#' @param keyword Character scalar used for naming the file
#'   (e.g., "species", "reef", "size", "dovs").
#' @param path Character scalar, directory where updated files should be saved.
#'   Default is `"data/lists/updates/"`.
#'
#' @return Invisibly returns the output filepath.
#'
#' @details
#' - Output filenames follow the pattern:
#'   - `"ltem_monitoring_<keyword>_YYYY-MM-DD.xlsx"` for species, reefs, sizes.
#'   - `"ltem_<keyword>_list_YYYY-MM-DD.xlsx"` if `<keyword>` already contains "list".
#' - The function creates the output directory if it does not exist.
#' - File format is always `.xlsx` using `writexl::write_xlsx()`.
#'
#' @examples
#' # Save updated species list
#' # save_ltem_list(species_list, "species")
#'
#' # Save updated reefs list
#' # save_ltem_list(reef_list, "reef")
#'
#' @keywords data-export LTEM save
## Save an updated LTEM list with an auto-generated date-stamped filename
save_ltem_list <- function(
    data,
    keyword,
    path = "data/lists/updates/"
) {
  # Ensure output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Today's date
  today <- format(Sys.Date(), "%Y-%m-%d")
  
  # Corrected case-insensitive detection
  if (str_detect(keyword, regex("list", ignore_case = TRUE))) {
    # e.g., size_list, dovs_reefs_list
    filename <- paste0("ltem_", keyword, "_", today, ".xlsx")
  } else {
    # General monitoring lists
    filename <- paste0("ltem_monitoring_", keyword, "_", today, ".xlsx")
  }
  
  out_path <- file.path(path, filename)
  
  writexl::write_xlsx(data, out_path)
  
  message("✔ File saved to: ", out_path)
  invisible(out_path)
}
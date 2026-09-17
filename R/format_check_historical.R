#' Map size IDs to canonical sizes (historical) and flag modifications
#'
#' Joins an LTEM historical dataset with a size reference to correct size
#' values and track original values and status.
#'
#' @param LTEM Historical LTEM dataframe with `Label`, `IDSize`, and `Size`.
#' @param IDSize Reference dataframe mapping (`Label`, `IDSize`) to `Size`.
#' @return Dataframe with updated `Size`, plus `Before` and `Status`.
#' @examples
#' # out <- sizeid(hist_ltem, id_size)
sizeid_hist <- function(LTEM, IDSize){
  return(sizeid(LTEM, IDSize))
}

#' Map reef names to IDs (historical) and flag modifications
#'
#' Joins LTEM historical data to a reef reference to correct `IDReef` and
#' track changes.
#'
#' @param LTEM Historical LTEM dataframe with `Reef`, `Habitat`, `IDReef`.
#' @param IDReef Reference dataframe with (`IDReef`, `Reef`, `Habitat`).
#' @return Dataframe with corrected `IDReef`, plus `IDBefore` and `Status`.
#' @examples
#' # out <- reefsid(hist_ltem, id_reef)
reefsid_hist <- function(LTEM, IDReef){
  return(reefsid(LTEM, IDReef))
}

#' Replace reef names with canonical names (historical)
#'
#' Joins by `IDReef` to bring reference `Reef`, comparing and updating
#' values while preserving previous value and status.
#'
#' @param LTEM Historical LTEM dataframe with `IDReef` and `Reef`.
#' @param IDReef Reference dataframe with (`IDReef`, `Reef`).
#' @return Dataframe with updated `Reef`, plus `ReefBefore` and `Status`.
#' @examples
#' # out <- reefsname(hist_ltem, id_reef)
reefsname_hist <- function(LTEM, IDReef){
  return(reefsname(LTEM, IDReef))
}

#' Build a flag table of modified records (historical)
#'
#' Filters rows whose `Status` contains "Modified_" and selects identifiers
#' and correction columns to report.
#'
#' @param LTEM_corrected Historical LTEM dataframe with correction columns.
#' @return Dataframe of flagged rows.
#' @examples
#' # flagged <- flags(hist_ltem_corrected)
flags_hist <- function(LTEM_corrected){
  return(flags(LTEM_corrected))
}

#' Extract unique transect entries (historical)
#'
#' De-duplicates transects by year and region with basic metadata.
#'
#' @param LTEM Historical LTEM dataframe.
#' @return Dataframe of unique transects.
#' @examples
#' # uniq <- unique_trnsct(hist_ltem)
unique_trnsct_hist <- function(LTEM){
  return(unique_trnsct(LTEM))
}

#' Compare INV and FISH transect coverage (historical)
#'
#' Merges invertebrate and fish transects to identify missing labels across
#' datasets.
#'
#' @param INV Historical invertebrate transects dataframe.
#' @param FISH Historical fish transects dataframe.
#' @return Dataframe with both labels per transect.
#' @examples
#' # cmp <- compare_trnsct(inv_hist, fish_hist)
compare_trnsct_hist <- function(INV, FISH) {
  return(compare_trnsct(INV, FISH))
}

#' Flag missing transects between datasets (historical)
#'
#' Marks transects missing in either INV or FISH and filters out those
#' present in both. Also writes an RDS snapshot to disk.
#'
#' @param COMPARED Output of `compare_trnsct()`.
#' @return Dataframe of missing transects grouped by region.
#' @examples
#' # miss <- flag_trnsct(cmp_hist)
flag_trnsct_hist <- function(COMPARED)  {
  m_transects <- flag_trnsct(COMPARED)
  readr::write_rds(m_transects, "data/missing_transects.RDS")
  return(m_transects)
}

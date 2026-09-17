#' Map Species IDs (historical) and flag modifications
#'
#' Joins historical LTEM data with a species catalog by `Label` and
#' `Species` to bring canonical `IDSpecies`, preserving additional columns
#' where available and tracking prior values and status.
#'
#' @param LTEM Historical LTEM dataframe with `Label`, `Species`, `IDSpecies`.
#' @param clean_sp_list Historical species catalog including trait columns.
#' @return Dataframe with corrected `IDSpecies`, plus `IDBefore` and `Status`.
#' @examples
#' # out <- speciesid(hist_ltem, species_list_hist)
speciesid_hist <- function(LTEM, clean_sp_list){
  return(speciesid(LTEM, clean_sp_list))
}



#' Replace species names by ID mapping (historical) and flag modifications
#'
#' Joins by `Label` and `IDSpecies` to bring canonical `Species` names from
#' the historical species catalog, tracking original value and modification
#' status; preserves associated trait columns.
#'
#' @param LTEM_IDs Historical LTEM dataframe with `Label`, `IDSpecies`.
#' @param clean_sp_list Historical species catalog.
#' @return Dataframe with updated `Species`, plus `SpeciesBefore` and `Status`.
#' @examples
#' # out <- speciesnames(hist_ids, species_list_hist)
speciesnames_hist <- function(LTEM_IDs, clean_sp_list){
  return(speciesnames(LTEM_IDs, clean_sp_list))
}
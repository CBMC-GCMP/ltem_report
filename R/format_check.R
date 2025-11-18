#' Map size IDs to canonical sizes and flag modifications
#'
#' Joins an LTEM dataset with a size reference table to replace mismatched
#' sizes, while tracking the original value and status.
#'
#' @param LTEM Data frame with at least `Label`, `IDSize`, and `Size`.
#' @param IDSize Data frame mapping (`Label`, `IDSize`) to canonical `Size`.
#' @return Data frame with updated `Size`, plus `Before` and `Status`.
#' @examples
#' # out <- sizeid(ltem, id_size)
sizeid <- function(LTEM,IDSize){
  ltem <- dplyr::left_join(LTEM, IDSize, by = c("Label","IDSize"), suffix = c(".x", ".y")) %>% 
    dplyr::rename(correct_size = Size.y,
                  Size = Size.x) %>% 
    dplyr::mutate(Before = Size,
                  Status = dplyr::if_else(Size == correct_size, "Correct", "Modified"),
                  Size = dplyr::if_else(Size == correct_size, Size, correct_size))
  return(ltem)
}


#' Map reef names to IDs and flag modifications
#'
#' Normalizes reef names (replacing spaces with underscores), then joins to
#' a reef reference to correct `IDReef` and track changes.
#'
#' @param LTEM Data frame with `Reef`, `Habitat`, and `IDReef`.
#' @param IDReef Data frame with (`IDReef`, `Reef`, `Habitat`).
#' @return Data frame with corrected `IDReef`, plus `IDBefore` and `Status`.
#' @examples
#' # out <- reefsid(ltem, id_reef)
reefsid <- function(LTEM,IDReef){
  
  LTEM <- LTEM %>% 
    dplyr::mutate(Reef = stringr::str_replace(Reef, " ", "_"))
  
  ltem <- dplyr::left_join(LTEM, IDReef[, c("IDReef", "Reef", "Habitat")], by = c("Reef", "Habitat"), suffix = c(".x", ".y")) %>% 
    dplyr::rename(correct_id = IDReef.y,
                  IDReef = IDReef.x) %>% 
    dplyr::mutate(IDBefore = IDReef,
                  Status = dplyr::case_when(IDReef == correct_id ~ "Correct",
                                             is.na(IDReef) ~ "Modified_ID",
                                             TRUE ~ "Modified_ID"),
                  IDReef = dplyr::if_else(is.na(IDReef), correct_id, IDReef),
                  IDReef = dplyr::if_else(IDReef == correct_id, IDReef, correct_id),
                  IDReef = dplyr::if_else(is.na(correct_id), IDBefore, IDReef))
  return(ltem)
}


#' Replace reef names with canonical names and flag modifications
#'
#' Joins by `IDReef` to bring the reference `Reef` name, comparing and
#' updating values while preserving the previous value and status.
#'
#' @param LTEM Data frame with `IDReef` and `Reef`.
#' @param IDReef Data frame with (`IDReef`, `Reef`).
#' @return Data frame with updated `Reef`, plus `ReefBefore` and `Status`.
#' @examples
#' # out <- reefsname(ltem, id_reef)
reefsname <- function(LTEM,IDReef){
  ltem <- dplyr::left_join(LTEM, IDReef[, c("IDReef", "Reef")], by = "IDReef", suffix = c(".x", ".y")) %>% 
    dplyr::rename(correct_reef = Reef.y,
                  Reef = Reef.x) %>% 
    dplyr::mutate(ReefBefore = Reef,
                  Status = dplyr::if_else(Reef == correct_reef, Status, "Modified_Reef"),
                  Reef = dplyr::if_else(Reef == correct_reef, Reef, correct_reef))
  return(ltem)
}


#' Extract unique transect entries
#'
#' Produces a de-duplicated list of transects by region/date/reef/depth and
#' observer.
#'
#' @param LTEM Data frame with transect metadata.
#' @return Data frame of unique transects ordered by date and location.
#' @examples
#' # uniq <- unique_trnsct(ltem)
unique_trnsct <- function(LTEM){
  transects <- LTEM %>% 
  dplyr::group_by(Region, Year, Month, Day, IDReef, Reef, Depth) %>% 
  dplyr::select(Label, Region, Year, Month, Day, IDReef, Reef, Depth, Transect, Observer) %>% 
  dplyr::arrange(Year, Month, Day, IDReef, Reef, Depth, Transect) %>% 
  unique()
  
  return(transects)
}


#' Compare INV and FISH transect coverage
#'
#' Merges invertebrate and fish transect tables to identify missing labels
#' per transect across both datasets.
#'
#' @param INV Data frame of INV transects.
#' @param FISH Data frame of FISH transects.
#' @return Data frame with both labels per transect.
#' @examples
#' # cmp <- compare_trnsct(inv, fish)
compare_trnsct <- function(INV, FISH) {
  keys <- c("Region","Year","Month","Day","IDReef","Reef","Habitat","Depth","Transect")
  join_keys <- intersect(keys, intersect(names(INV), names(FISH)))
  t_merge <- dplyr::full_join(INV, FISH, by = join_keys, suffix = c(".inv",".pec")) %>% 
    dplyr::mutate(
      Label_inv = dplyr::coalesce(rlang::.data[["Label.inv"]], rlang::.data[["Label.x"]], rlang::.data[["Label"]]),
      Label_pec = dplyr::coalesce(rlang::.data[["Label.pec"]], rlang::.data[["Label.y"]], rlang::.data[["Label"]])
    ) %>%
    dplyr::select(-dplyr::any_of(c("Label.inv","Label.pec","Label.x","Label.y","Label")))
  
  return(t_merge)
}


#' Flag missing transects between datasets
#'
#' Marks transects missing in either INV or FISH and filters out those
#' present in both.
#'
#' @param COMPARED Output of `compare_trnsct()`.
#' @return Data frame of missing transects grouped by region.
#' @examples
#' # miss <- flag_trnsct(cmp)
flag_trnsct <- function(COMPARED)  {
  m_transects <- COMPARED %>% 
  dplyr::mutate(Flag = dplyr::if_else(is.na(Label_inv), "missing_INV", "CORRECT"),
                Flag = dplyr::if_else(is.na(Label_pec), "missing_PEC", Flag)
  ) %>% 
  dplyr::group_by(Region) %>% 
  dplyr::filter(Flag != "CORRECT")
  # write_rds(COMPARED,"data/missing_transects.RDS") 
  return(m_transects)
}

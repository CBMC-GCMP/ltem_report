#' Map Species IDs and flag modifications
#'
#' Joins LTEM data with a species catalog by `Label` and `Species` to bring
#' canonical `IDSpecies`, tracking the prior value and a modification status.
#'
#' @param LTEM Data frame with columns `Label`, `Species`, and `IDSpecies`.
#' @param clean_sp_list Species catalog with `Label`, `IDSpecies`, `Species`.
#' @return Data frame with corrected `IDSpecies`, plus `IDBefore` and `Status`.
#' @examples
#' # out <- speciesid(ltem, species_list)
speciesid <- function(LTEM, clean_sp_list){

        merge_df <- dplyr::left_join(
          LTEM,
          dplyr::select(clean_sp_list, Label, IDSpecies, Species, dplyr::any_of(c("A_ord","B_pen","TrophicLevel","Functional_groups"))),
          by = c("Label", "Species"),
          suffix = c(".x", ".y")
        ) %>%
          dplyr::rename(correct_id = IDSpecies.y,
                        IDSpecies = IDSpecies.x) %>%
          dplyr::mutate(IDBefore = IDSpecies,
                         Status = dplyr::if_else(IDSpecies == correct_id, "Correct", "Modified_IDSp"),
                         IDSpecies = dplyr::if_else(is.na(IDSpecies), correct_id, IDSpecies),
                         IDSpecies = dplyr::if_else(IDSpecies == correct_id, IDSpecies, correct_id),
                         IDSpecies = dplyr::if_else(is.na(correct_id), IDBefore, IDSpecies))

        return(merge_df)

    }


    #' Replace species names by ID mapping and flag modifications
    #'
    #' Joins by `Label` and `IDSpecies` to bring canonical `Species` names from
    #' the species catalog, tracking original value and a modification status.
    #'
    #' @param LTEM_IDs Data frame with `Label` and `IDSpecies` columns.
    #' @param clean_sp_list Species catalog with `Label`, `IDSpecies`, `Species`.
    #' @return Data frame with updated `Species`, plus `SpeciesBefore` and `Status`.
    #' @examples
    #' # out <- speciesnames(ltem_ids, species_list)
    speciesnames <- function(LTEM_IDs, clean_sp_list){
      
    
      LTEM <- LTEM_IDs
    
        
        merge_df <- dplyr::left_join(
          LTEM,
          dplyr::select(clean_sp_list, Label, IDSpecies, Species, dplyr::any_of(c("A_ord","B_pen","TrophicLevel","Functional_groups"))),
          by = c("Label", "IDSpecies"),
          suffix = c(".x", ".y")
        ) %>% 
          dplyr::rename(correct_sp = Species.y,
                        Species = Species.x) %>%
          dplyr::mutate(SpeciesBefore = Species,
                         Status = dplyr::if_else(Species == correct_sp, "Correct", "Modified_sp"),
                         Species = dplyr::if_else(Species == correct_sp, Species, correct_sp),
                         Species = dplyr::if_else(is.na(correct_sp), SpeciesBefore, Species))

        return(merge_df)
      
    }

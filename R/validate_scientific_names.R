#' Validate and reconcile scientific names with external sources
    #'
    #' Merges resolved names into a cleaned species list and, depending on the
    #' label, validates fish names via rfishbase or resolves invertebrate names
    #' via WoRMS IDs. Also writes flag reports to `data/flag_reports/` and
    #' species list updates to `data/lists/updates/` when applicable.
    #'
    #' @param clean_md Data frame of cleaned species metadata, including `Species`.
    #' @param resolved_names Data frame with columns `Species` and
    #'   `resolved_scientific_name` to apply.
    #' @param species_list Species catalog to update.
    #' @param Label Character scalar, either "fish" or "inv".
    #' @return For `Label = "fish"`, returns the updated `species_list` with
    #'   validated names. For `Label = "inv"`, returns a data frame of
    #'   `valid_names` containing `IDSpecies` and `Species` (resolved names).
    #' @examples
    #' # out <- clean_validation(clean_md, resolved, species_list, "fish")
    clean_validation <- function(clean_md, resolved_names, species_list, Label) {
      
      merge_md <-
        merge(clean_md, resolved_names[, c("Species", "resolved_scientific_name")],
              by = "Species", all.x = TRUE) %>%
        mutate(resolved_scientific_name= ifelse(is.na(resolved_scientific_name), 
                                                Species, resolved_scientific_name),
               SpeciesBefore = Species,
               Status = ifelse(
                 Species == resolved_scientific_name,
                 "correct_sp",
                 "Modified_sp"
               ),
               Species = ifelse(
                 Species == resolved_scientific_name,
                 Species,
                 resolved_scientific_name
               )
        )
      
      dir.create("data/flag_reports", showWarnings = F)
      dir.create("data/lists/updates", showWarnings = F)
      
      if(Label=="fish"){
        valid_names <- merge_md %>%
          mutate(
            New =  rfishbase::validate_names(Species),
            correct_sp = ifelse(Species == New, Species, New),
            Status = ifelse(Species == correct_sp, Status, "Modified_sp")
          ) %>%
          select(-New)
        
        
        flags <- valid_names %>%
          filter(str_detect(Status, "Modified_")) %>% 
          select(IDSpecies, Species, SpeciesBefore, correct_sp, Status) 
        writexl::write_xlsx(flags, "data/flag_reports/fish_validated_sci-names.xlsx")
        
        valid_names <- valid_names %>%
          select(IDSpecies, correct_sp) 
        
        species_list <- merge(species_list, valid_names, by="IDSpecies", all=T) %>% 
          mutate(correct_sp= ifelse(is.na(correct_sp),Species, correct_sp ),
                 Species= ifelse(Species==correct_sp, Species, correct_sp)) %>% 
          select(-c(correct_sp))
        
        write.csv(species_list, "data/lists/updates/fish_updated_names.csv")
        species_list <- species_list 
        return(species_list)
        
      }else if(Label=="inv"){
        flags <- merge_md  %>%
          filter(str_detect(Status, "Modified_")) %>% 
          select(IDSpecies, Species, SpeciesBefore, resolved_scientific_name, Status) 
        writexl::write_xlsx(flags, "data/flag_reports/inv_resolved_sci-names.xlsx")
        
        valid_names <- merge_md  %>%
          select(IDSpecies, Species) 
        return(valid_names)
      }else{
        stop("Please specify 'fish' or 'inv'")
      }
      
    }

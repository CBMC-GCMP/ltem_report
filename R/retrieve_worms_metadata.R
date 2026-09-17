#Retrieve worms IDs using taxize:: get_wormsid
## INVs Only
#' Retrieve WoRMS IDs for invertebrate species
#'
#' Uses `taxize::get_wormsid()` to retrieve WoRMS AphiaIDs for species names.
#'
#' @param correct_sp Data frame with columns `IDSpecies` and `Species`.
#' @return Data frame of WoRMS IDs as returned by `taxize::get_wormsid()`.
#' @examples
#' # ids <- check_worms(correct_sp)
check_worms <- function(correct_sp){
  spp_inv <- data.frame(IDSpecies= correct_sp$IDSpecies,
                        Species= stringr::str_replace(correct_sp$Species,
                                                      " ", "+"))
  wormsID <- taxize:: get_wormsid(spp_inv$Species)
  wormsID <- data.frame(wormsID)
  return(wormsID)
}

# Transform WoRMS IDs into validated scientific names, and return them
# to our PLoS

#' Convert WoRMS IDs to validated names and update species list
#'
#' Given AphiaIDs, fetch records from WoRMS and reconcile against the
#' existing species list to update canonical names.
#'
#' @param wormsID Data frame (from `check_worms()`) containing an `ids` column.
#' @param clean_md Data frame of species with `IDSpecies` and `Species`.
#' @param species_list Species catalog to update.
#' @return Updated `species_list` data frame with validated names written to
#'   `data/lists/updates/`.
#' @examples
#' # updated <- worms_format(worms_ids, clean_md, species_list)
worms_format <- function(wormsID, clean_md, species_list){
  
  wormsID = as.numeric(wormsID$ids)
  
  worms_names <- worms::wormsbyid(wormsID,
                                  verbose = FALSE,
                                  ids = TRUE,
                                  sleep_btw_chunks_in_sec = 0.01)
  
  ## We clean our worms_names df, and add our IDSpecies
  
  clean_sp_list <- clean_md %>% 
    mutate(Species= str_replace_all(Species, "[+]", " "))
  
  
  
  valid_names<- worms_names %>% 
    mutate(IDSpecies= ifelse(clean_sp_list$Species== scientificname,
                             clean_sp_list$IDSpecies, "")) %>% 
    select(IDSpecies, valid_name, valid_AphiaID)
  
  
  
  ## We merge our valid_names with our clean PLoS
  
  species_list <- merge(species_list, valid_names, by="IDSpecies", all=T) %>% 
    mutate(valid_name= ifelse(is.na(valid_name),Species, valid_name ),
           Species= ifelse(Species==valid_name, Species, valid_name)) %>% 
    select(-c(valid_name))
  writexl::write_xlsx(species_list, paste0("data/lists/updates/ltem_monitoring_species_", Sys.Date(),".xlsx"))
  species_list <- species_list  
  return(species_list)
  
}

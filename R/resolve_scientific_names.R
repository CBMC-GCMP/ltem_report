##It is easier to work with a peripherical list of Species (PLoS)

### First step is to correct the spelling

##We load WoRMS and Fishbase as our sources

##Then, we apply a taxsize function: gna_verifier, for matching correct
## spellings of Scientific Names, in three parts:


##Both INVs & FISH

#' Resolve scientific names against WoRMS and FishBase
#'
#' Uses `taxize::gna_verifier()` to retrieve best matches for scientific
#' names from WoRMS (and FishBase for fish). Filters out genus-only matches
#' and returns a two-column data frame with input `Species` and a
#' `resolved_scientific_name`.
#'
#' @param clean_spp Data frame with a `Species` column of candidate names.
#' @param Label Character scalar, either "fish" or "inv". Determines data
#'   sources and filtering.
#' @return Data frame with columns `Species` and `resolved_scientific_name`.
#' @examples
#' # rn <- resolve_names(clean_spp, "fish")
resolve_names <- function(clean_spp, Label){
  
  #FISH
  if(Label== "fish"){
    
    ## Part I: 
    #Retrieve best species matches from WoRMS and Fishbase      
    sources <- c(worms = 9, fishbase = 155)
    resolved_names <- sources %>% 
      map(~ taxize::gna_verifier(data_source_ids = .x, 
                                names = clean_spp$Species,
                                best_match_only = T, 
                                fields = c("all"), 
                                canonical = T ))
    
    ## Part II:
    #Leave only the correct results, and format it as data.frame
    # Check the structure of the output and adapt accordingly
    
    resolved_names_df <- resolved_names %>% 
      map(~ {
        # Check if the output has the expected columns
        if("match_value" %in% names(.x)) {
          # Old format
          .x %>% 
            filter(!match_value %in% c("Could only match genus") & 
                     str_count(matched_name2, "\\w+") >= 2) %>% 
            select(user_supplied_name, matched_name2, taxon_id)
        } else {
          # New format - based on the actual column names we see
          # Filter out genus-only matches and ensure we have at least genus and species
          .x %>% 
            filter(!is.na(matchedCanonicalSimple) & 
                     str_count(matchedCanonicalSimple, "\\w+") >= 2) %>% 
            select(user_supplied_name = submittedName, 
                   matched_name2 = matchedCanonicalSimple,
                   taxon_id = matchedNameID)
        }
      }) %>% 
      reduce(full_join, by = "user_supplied_name") %>% 
      set_names(c("user_supplied_name", "worms_sci_name",
                  "worms_id", "fishbase_sci_name", "fishbase_id"))
    
    
    ## Part III:
    #Filter absent values (NAs)
    
    resolved_names <- resolved_names_df %>% 
      mutate(resolved_scientific_name = case_when(
        !is.na(worms_sci_name) ~ worms_sci_name,
        !is.na(fishbase_sci_name) ~ fishbase_sci_name,
        TRUE ~ user_supplied_name
      )) %>% 
      select(user_supplied_name, 
             resolved_scientific_name, 
             everything())%>% 
      rename(Species= user_supplied_name) %>% 
      select(Species, resolved_scientific_name)
    
    
    #INVs
  }else if(Label=="inv"){
    ##Part I:
    sources <- c(worms = 9)
    resolved_names <- sources %>% 
      map(~ taxize::gna_verifier(data_source_ids = .x, 
                                names = clean_spp$Species,
                                best_match_only = T, 
                                fields = c("all"), 
                                canonical = T ))
    
    #Part II:
    resolved_names_df <- resolved_names %>% 
      map(~ {
        # Check if the output has the expected columns
        if("match_value" %in% names(.x)) {
          # Old format
          .x %>% 
            filter(!match_value %in% c("Could only match genus") &
                     str_count(matched_name2, "\\w+") >= 2) %>% 
            select(user_supplied_name, matched_name2, taxon_id)
        } else {
          # New format - based on the actual column names we see
          # Filter out genus-only matches and ensure we have at least genus and species
          .x %>% 
            filter(!is.na(matchedCanonicalSimple) & 
                     str_count(matchedCanonicalSimple, "\\w+") >= 2) %>% 
            select(user_supplied_name = submittedName, 
                   matched_name2 = matchedCanonicalSimple,
                   taxon_id = matchedNameID)
        }
      }) %>% 
      reduce(full_join, by = "user_supplied_name") %>% 
      set_names(c("user_supplied_name",
                  "worms_sci_name",
                  "worms_id"))
    
    #Part III:
    resolved_names <- resolved_names_df %>% 
      mutate(resolved_scientific_name = case_when(
        !is.na(worms_sci_name) ~ worms_sci_name,
        TRUE ~ user_supplied_name
      )) %>% 
      select(user_supplied_name, 
             resolved_scientific_name, 
             everything())%>% 
      rename(Species= user_supplied_name) %>% 
      select(Species, resolved_scientific_name)
  } else{
    stop("Please specify 'fish' or 'inv'")
  }
  
  return(resolved_names)
}

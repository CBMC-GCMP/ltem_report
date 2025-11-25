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
    input_names <- unique(stringr::str_squish(clean_spp$Species))
    input_names <- input_names[!is.na(input_names) & nzchar(input_names)]
    chunks <- split(input_names, ceiling(seq_along(input_names) / 50))
    fetch_source <- function(src_id){
      purrr::map_dfr(chunks, function(ch){
        res <- purrr::possibly(taxize::gna_verifier, otherwise = NULL)(
          data_source_ids = src_id,
          names = ch,
          best_match_only = TRUE,
          fields = c("all"),
          canonical = TRUE
        )
        if (is.null(res)) {
          purrr::map_dfr(ch, function(nm){
            res1 <- purrr::possibly(taxize::gna_verifier, otherwise = NULL)(
              data_source_ids = src_id,
              names = nm,
              best_match_only = TRUE,
              fields = c("all"),
              canonical = TRUE
            )
            if (is.null(res1)) dplyr::tibble() else res1
          })
        } else {
          res
        }
      })
    }
    resolved_names <- purrr::map(sources, fetch_source)
    
    ## Part II:
    #Leave only the correct results, and format it as data.frame
    # Check the structure of the output and adapt accordingly
    
    resolved_names_df <- resolved_names %>% 
      map(~ {
        if (is.null(.x) || (is.data.frame(.x) && nrow(.x) == 0)) {
          return(dplyr::tibble(user_supplied_name = character(), matched_name2 = character(), taxon_id = character()))
        }
        # Check if the output has the expected columns
        if("match_value" %in% names(.x)) {
          # Old format
          .x %>% 
            filter(!match_value %in% c("Could only match genus") & 
                     str_count(matched_name2, "\\w+") >= 2) %>% 
            select(user_supplied_name, matched_name2, taxon_id)
        } else if (all(c("submittedName","matchedCanonicalSimple","matchedNameID") %in% names(.x))) {
          # New format - based on the actual column names we see
          # Filter out genus-only matches and ensure we have at least genus and species
          .x %>% 
            filter(!is.na(matchedCanonicalSimple) & 
                     str_count(matchedCanonicalSimple, "\\w+") >= 2) %>% 
            select(user_supplied_name = submittedName, 
                   matched_name2 = matchedCanonicalSimple,
                   taxon_id = matchedNameID)
        } else {
          dplyr::tibble(user_supplied_name = character(), matched_name2 = character(), taxon_id = character())
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
    input_names <- unique(stringr::str_squish(clean_spp$Species))
    input_names <- input_names[!is.na(input_names) & nzchar(input_names)]
    chunks <- split(input_names, ceiling(seq_along(input_names) / 50))
    fetch_source <- function(src_id){
      purrr::map_dfr(chunks, function(ch){
        res <- purrr::possibly(taxize::gna_verifier, otherwise = NULL)(
          data_source_ids = src_id,
          names = ch,
          best_match_only = TRUE,
          fields = c("all"),
          canonical = TRUE
        )
        if (is.null(res)) {
          purrr::map_dfr(ch, function(nm){
            res1 <- purrr::possibly(taxize::gna_verifier, otherwise = NULL)(
              data_source_ids = src_id,
              names = nm,
              best_match_only = TRUE,
              fields = c("all"),
              canonical = TRUE
            )
            if (is.null(res1)) dplyr::tibble() else res1
          })
        } else {
          res
        }
      })
    }
    resolved_names <- purrr::map(sources, fetch_source)
    
    #Part II:
    resolved_names_df <- resolved_names %>% 
      map(~ {
        if (is.null(.x) || (is.data.frame(.x) && nrow(.x) == 0)) {
          return(dplyr::tibble(user_supplied_name = character(), matched_name2 = character(), taxon_id = character()))
        }
        # Check if the output has the expected columns
        if("match_value" %in% names(.x)) {
          # Old format
          .x %>% 
            filter(!match_value %in% c("Could only match genus") &
                     str_count(matched_name2, "\\w+") >= 2) %>% 
            select(user_supplied_name, matched_name2, taxon_id)
        } else if (all(c("submittedName","matchedCanonicalSimple","matchedNameID") %in% names(.x))) {
          # New format - based on the actual column names we see
          # Filter out genus-only matches and ensure we have at least genus and species
          .x %>% 
            filter(!is.na(matchedCanonicalSimple) & 
                     str_count(matchedCanonicalSimple, "\\w+") >= 2) %>% 
            select(user_supplied_name = submittedName, 
                   matched_name2 = matchedCanonicalSimple,
                   taxon_id = matchedNameID)
        } else {
          dplyr::tibble(user_supplied_name = character(), matched_name2 = character(), taxon_id = character())
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

### Apply filters for removing all species entries that contain "sp" or "spp" 
### Including Species with just genus or families
##Filtering just Fish data

##Both INVs & FISH
#' Clean species names by filtering ambiguous entries
#'
#' Applies filters to keep only binomial species names and remove ambiguous
#' entries such as "sp", "spp", and "cf." tokens. Behavior differs for fish
#' vs. invertebrates based on the `Label` column in the input.
#'
#' @param species_list A data.frame/tibble with at least columns `Label` and
#'   `Species` (and optionally `IDSpecies`).
#' @param Label Character scalar, either "fish" or "inv". Selects records by
#'   `Label` ("PEC" for fish, "INV" for invertebrates) and chooses the filter set.
#' @return A data.frame/tibble with filtered rows. For `Label = "fish"`, returns
#'   only `Label`, `IDSpecies`, and `Species`. For `Label = "inv"`, preserves
#'   input columns.
#' @details
#' - Fish ("PEC"): keeps names containing a space and excludes " sp" or "spp".
#' - Invertebrates ("INV"): keeps names containing a space and excludes
#'   " sp1"..." sp9", "spp", and " cf"/" cf.".
#' @examples
#' # cleaned <- clean_spp(species_list, "fish")
#' # cleaned_inv <- clean_spp(species_list, "inv")
#' @keywords data-cleaning
clean_spp <- function(species_list, Label){
  if(Label=="fish"){
    clean_md <-  species_list %>% 
      select(Label, IDSpecies, Species ) %>% 
      filter (Label== "PEC") %>%
      filter(str_detect(Species, " " )) %>% 
      filter(!str_detect(Species, " sp | spp" ))
  } else if(Label=="inv"){
    clean_md<-  species_list %>% 
      filter (Label== "INV") %>% 
      filter(str_detect(Species, " " )) %>%
      filter(!str_detect(Species, " sp1| sp2| sp3| sp4| sp5| sp6| sp7| sp8| sp9| spp" )) %>% 
      filter(!str_detect(Species, " cf| cf."))    
  }else {
    stop("Please specify 'fish' or 'inv'")
  }
  
  return(clean_md)
}

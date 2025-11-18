library(tidyverse)
library(readxl)
library(namechecker)


# FISH --------------------------------------------------------------------
## Load custom functions

lapply(list.files(path="R/", pattern = ".R", full.names = T), source)

pipeline_taxonomy_update <- function(cfg) {
  lists_dir <- cfg$lists_dir
  outputs_dir <- cfg$outputs_dir
  if (is.null(lists_dir) || is.null(outputs_dir)) stop("cfg must include 'lists_dir' and 'outputs_dir'")
  review_dir <- file.path(outputs_dir, "review", "taxonomy")
  dir.create(review_dir, showWarnings = FALSE, recursive = TRUE)
  spp_files <- list.files(lists_dir, pattern = "^ltem_monitoring_species.*\\.xlsx$", full.names = TRUE)
  if (length(spp_files) == 0) stop("No species list files found in: ", lists_dir)
  spp_path <- spp_files[which.max(file.info(spp_files)$mtime)]
  species_list <- readxl::read_xlsx(spp_path)
  original <- species_list[, c("IDSpecies","Species")]
  if (requireNamespace("namechecker", quietly = TRUE)) {
    valid <- namechecker::get_valid_names(species_list$Species)
    valid <- valid[, c("input_name","valid_name")]
    valid$valid_name[is.na(valid$valid_name)] <- valid$input_name
    tmp <- merge(species_list, valid, by.x = "Species", by.y = "input_name", all.x = TRUE)
    idx_map <- !is.na(tmp$valid_name)
    tmp$Species[idx_map] <- tmp$valid_name[idx_map]
    tmp$valid_name <- NULL
    species_list <- tmp
  }
  fish_updated <- species_list
  fish_clean <- tryCatch(clean_spp(species_list, "fish"), error = function(e) species_list[0, , drop = FALSE])
  resolved_fish <- tryCatch(resolve_names(fish_clean, "fish"), error = function(e) NULL)
  if (!is.null(resolved_fish) && nrow(fish_clean) > 0) {
    fish_updated <- tryCatch(clean_validation(fish_clean, resolved_fish, species_list, "fish"), error = function(e) species_list)
  }
  inv_clean <- tryCatch(clean_spp(fish_updated, "inv"), error = function(e) fish_updated[0, , drop = FALSE])
  resolved_inv <- tryCatch(resolve_names(inv_clean, "inv"), error = function(e) NULL)
  inv_valid <- NULL
  if (!is.null(resolved_inv) && nrow(inv_clean) > 0) {
    inv_valid <- tryCatch(clean_validation(inv_clean, resolved_inv, fish_updated, "inv"), error = function(e) NULL)
  }
  final_species <- fish_updated
  if (!is.null(inv_valid) && is.data.frame(inv_valid) && nrow(inv_valid) > 0) {
    wormsID <- tryCatch(check_worms(inv_valid), error = function(e) NULL)
    if (!is.null(wormsID)) {
      final_species <- tryCatch(worms_format(wormsID, inv_clean, fish_updated), error = function(e) fish_updated)
    }
  }
  flags_path <- NA_character_
  if (all(c("IDSpecies","Species") %in% names(final_species))) {
    before <- original
    after <- final_species[, c("IDSpecies","Species")]
    m <- merge(before, after, by = "IDSpecies", suffixes = c("Before",""))
    m$Status <- ifelse(is.na(m$SpeciesBefore) | m$SpeciesBefore == m$Species, "correct_sp", "Modified_sp")
    flags <- m[, c("IDSpecies","Species","SpeciesBefore","Status")]
    if (exists("save_stage", mode = "function")) {
      flags_path <- save_stage(flags, review_dir, 0, "taxonomy-flags", fmt = "csv")
    } else {
      flags_path <- file.path(review_dir, sprintf("stage-00_taxonomy-flags_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
      utils::write.csv(flags, flags_path, row.names = FALSE)
    }
  }
  out_catalog <- file.path(lists_dir, sprintf("ltem_monitoring_species_%s.xlsx", Sys.Date()))
  if (requireNamespace("writexl", quietly = TRUE)) {
    writexl::write_xlsx(final_species, out_catalog)
  } else {
    out_catalog <- sub("\\.xlsx$", ".csv", out_catalog)
    utils::write.csv(final_species, out_catalog, row.names = FALSE)
  }
  if (exists("stage_log_en_es", mode = "function")) {
    rows_updated <- tryCatch({
      af <- merge(original, final_species[, c("IDSpecies","Species")], by = "IDSpecies", all.x = TRUE)
      sum(af$Species.x != af$Species.y, na.rm = TRUE)
    }, error = function(e) NA_integer_)
    stage_log_en_es(0,
                    "Taxonomy update applied (namechecker + FishBase/WoRMS)",
                    "Actualización de taxonomía aplicada (namechecker + FishBase/WoRMS)",
                    review_dir,
                    rows_updated,
                    integer(0),
                    "Proceed to Stage 1.",
                    "Procede a la Etapa 1.")
  }
  invisible(list(species_catalog_path = out_catalog,
                 species_catalog_source = spp_path,
                 flags_path = flags_path,
                 species_list = final_species))
}

if (FALSE) {
#Peripheral list of species (PLoS)
species_list <- read_xlsx("data/lists/updates/ltem_monitoring_species_2024-04-23.xlsx")



valid.names <- get_valid_names(species_list$Species) |> 
  select(input_name, valid_name)|> 
  mutate(valid_name=ifelse(is.na(valid_name), input_name, valid_name))

species_list <- merge(species_list, valid.names, by.x = "Species", by.y = "input_name", all.x = T) |> 
                select(-Species) |> 
                rename(Species=valid_name) |> 
  select(IDSpecies, Species, everything())


writexl::write_xlsx(species_list, "data/lists/updates/ltem_monitoring_species_2025_04-21.xlsx")

### Apply filters for removing all species entries that contain "sp" or "spp" 
### Including Species with just genus or families
##Filtering just Fish data

#Function used: clean_spp()

clean_md <-clean_spp (species_list, "fish")

# Scientific Names Correction 

# Connect to WoRMS and fishbase API using complentary functions
# Function resolve_names() searches for species 
# scientific names correct spelling

resolved_names <- resolve_names(clean_md, "fish")



# Names Validation 

# Validates current scientific names and automatically updates them in our PLoS

# Function used: clean_validation()

fish_validated <- clean_validation(clean_md, resolved_names, species_list, "fish")

rm(clean_md, resolved_names)

writexl::write_xlsx(fish_validated, "data/lists/updates/ltem_monitoring_species_2025-04-21.xlsx")
# INVERTEBRATES -----------------------------------------------------------

#Peripherical List of Species (PLoS), with fish sci-names corrected

clean_spp <- clean_spp(fish_validated, "inv")


#Check possible misspellings in PLoS
resolved_names <- resolve_names(clean_spp, "inv")

#Merge with PLoS replacing sci-names, and generating flags
clean_sp_list <- clean_validation(clean_spp, resolved_names, Label="inv")





# WoRMS validation 

#Exploratory checkup of sci-names, in case manual corrections necessary 

wormsID <- check_worms(clean_sp_list)

## If any scientific is still unaccepted, we can manually check its status

#Console: If manual input is required, always select the species entry
#with an Accepted Status



## All species that display INVALID status, require manual correction.
## You can do so by replacing the old scientific name string, with a new one:
## For example:
clean_sp_list<- clean_sp_list %>%
  mutate(Species = recode(Species, "Aplysina aztecus" = "Aplysina azteca",
  "Echinaster tenuispina"="Echinaster (Othilia) tenuispina",
  "Mycale ramulosa"="Mycale (Zygomycale) ramulosa",
  "Pacifigorgia englemanni"="Pacifigorgia englemanii",
  "Thais planospira"="Tribulus planospira" ))


# The replacements above may vary, always check for new invalid species



## If all species displayed a VALID status, we then retrieve updated sci-names
## from WoRMS, and replace them in our PLoS

ltem_species <- worms_format(wormsID, clean_spp, fish_validated)
view(ltem_species)
rm(clean_sp_list,clean_spp,fish_validated,resolved_names, species_list, wormsID)

writexl::write_xlsx(ltem_species, "data/lists/updates/ltem_monitoring_species_2025-04-21.xlsx")
}

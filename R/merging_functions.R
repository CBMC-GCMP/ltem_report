#Merge monitoring database with species catalog

#' Merge database with species catalog by species identifiers
#'
#' @param db Data frame of monitoring records including `IDSpecies`, `Species`, `Label`.
#' @param spp Species catalog with matching keys.
#' @return Merged data frame.
mg_spp <-function(db, spp){
db <- merge(db, spp, by = c("IDSpecies", "Species", "Label"))
return(db)
}

#merge 

#' Merge database with reef reference by region/reef identifiers
#'
#' @param db Data frame of monitoring records including `Region`, `IDReef`, `Reef`.
#' @param reefs Reef reference table with matching keys.
#' @return Merged data frame.
mg_reef <-function(db, reefs){
db <- merge(db, reefs, by = c("Region", "IDReef", "Reef"))
return(db)
}

#Fish Biomass

#' Compute fish biomass and trophic level factor
#'
#' @param ltem Data frame of fish observations with `Quantity`, `A_ord`, `Size`, `B_pen`, and `Area`.
#' @return Data frame with added `Biomass` and `TrophicLevelF` columns.
biomass <-function(ltem){
ltem_biomass <- ltem %>% 
  mutate( Biomass = (Quantity * A_ord * (Size^B_pen))/(Area * 100),
          TrophicLevelF = cut(as.numeric(TrophicLevel), 
                              breaks = c(2, 2.5, 3, 3.5, 4, 4.6), 
                              labels = c("2-2.5", "2.5-3", "3-3.5", "3.5-4", "4-4.5"), 
                              right = FALSE)
  )
return(ltem_biomass)
}

#Invertebrate Abundance

#' Compute top invertebrate abundance per MPA
#'
#' @param db_sf Data frame of invertebrate observations with `Label`, `IMPA`/`MPA`, `Reef`, `Transect`, `Species`, `Quantity`.
#' @return Data frame of top species by abundance per MPA, with `Species` reordered within MPA.
abundance <-function(db_sf){
ltem_abundance <- db_sf %>% 
  filter(Label == "INV") %>% 
  group_by(IMPA, Reef, Transect, Species) %>%
  summarise(Abundance = sum(Quantity)) %>%
  group_by(MPA,Species) %>% 
  summarise(Abundance = mean(Abundance)) %>%
  group_by(MPA) %>%
  top_n(10, Abundance) %>%
  ungroup() %>% 
  mutate(Species= reorder_within(Species, Abundance, MPA))

return(ltem_abundance)
}

library(googlesheets4)
library(tidyverse)

#Custom functions
lapply(list.files(path="R/", pattern = ".R", full.names = T), source)

#Connect to AWS server where the db is stored
source("connect_db/dbconnect-source.R")

dbListTables(ltem_db) #Generates a list of tables

#Load LTEM historical db
historical <- tbl(ltem_db,
              "ltem_historical_database") %>% 
  collect()

#Load cleaned new data for LTEM

db <- readRDS("outputs/temp/ltem_clean_2025-11-21.RDS") %>% 
  mutate(bleaching_coverage=NA,
         IDSpecies=ifelse(Species=="Sargocentron suborbitale", 364, IDSpecies)) %>% 
  select(Label, everything())


#Load latest files
species <- load_latest_ltem_list("species") %>% 
  mutate(IDSpecies=as.numeric(IDSpecies))

reefs <- load_latest_ltem_list("reefs")%>% 
  mutate(IDReef=as.numeric(IDReef))




# Merge with LTEM db

ltem_sp <- left_join(db %>% select(-c(Label,Species)), species) 

ltem_meta <- merge(ltem_sp, reefs, by=c("IDReef", "Reef", "Region"), all.x = T) %>% 
  select(-Habitat.y) %>% 
  rename(Habitat=Habitat.x)

#Add biomass and standarize columns

ltem_new<- ltem_meta %>% 
  mutate(A_ord = as.numeric(A_ord), 
         B_pen = as.numeric(B_pen), 
         Size=as.numeric(Size),
         Biomass = (Quantity * A_ord * (Size^B_pen))/(Area * 100)) %>% 
  mutate(TrophicGroup = factor(TrophicGroup, 
                               levels = c("Piscivoro", 
                                          "Carnivoro", 
                                          "Herbivoro", 
                                          "Zooplanctivoro")), 
         Region = factor(Region),
         TrophicLevelF = cut(as.numeric(TrophicLevel), 
                             breaks = c(2, 2.5, 3, 3.5, 4, 4.6), 
                             labels = c("2-2.5", "2.5-3", "3-3.5", "3.5-4", "4-4.5"), 
                             right = FALSE)
  ) %>% 
  mutate(Depth2= case_when( Depth <= 10 ~ "Shallow",
                            Depth >= 15 ~ "Deep"),
         Degree= round(Latitude,0),
         Genus=Species) %>%
  separate(Genus, c("Genus", NA), sep=" ", fill="left")%>% 
  
  # select( -Observer) %>% 
  arrange(Label, Year, Region, Reef, Transect, Depth, Species) %>% 
  select(any_of(names(historical)), everything(), -c(IDSize, valid_AphiaID))

dir.create("outputs/cleaned/", showWarnings = F)

#Save current LTEM db
saveRDS(ltem_new, "outputs/cleaned/ltem_2025-11-21.RDS")



# Update Historical db ----------------------------------------------------

#Load historical db



# Add new entries
updated <- rbind(historical %>% select(-row_names), ltem_new) 


#Save updated LTEM db
saveRDS(updated, "outputs/historical/ltem_historic_updated_2025-11-20.RDS")


library(tidyverse)
library(readxl)

#Load raw LTEM db's with missing columns or data


raw1 <- read_xlsx("data/raw/2025/loreto_29-31_oct/Ltem_loreto_ARI.xlsx") %>% 
  select(-Size)

raw2 <- read_excel("data/raw/2025/loreto_29-31_oct/2511_CaboPulmo - 20251107 - DB_REGISTRO_MONITOREO_2025_PFA - Alex.xlsx") %>% 
  # glimpse() %>% 
  select(any_of(names(ari))) %>% 
  select(-c(IDReef, IDSpecies)) %>% 
  # glimpse() %>% 
  mutate(Species=recode(Species, 
                        "Myxilla incrustans"="Myxilla (Myxilla) incrustans",
                        "Macrorhynchia nuttingi"="Aglaophenia whiteleggei",
                        "Lobatus galeatus"="Titanostrombus galeatus",
                        "Hyotissa solida"="Hyotissa hyotis"))

#Load size IDs db

sizes <- read.csv("data/lists/ltem_size_list.csv")

#Load Monitoring Reefs db
reefs <- read_excel("data/lists/ltem_monitoring_reefs_2025-05-06.xlsx") %>% 
  select(IDReef, Region, Reef, Habitat)

#Load Monitoring Species db
species <- read_excel("data/lists/ltem_monitoring_species_2025_04-21.xlsx") %>% 
  select(IDSpecies, Species)
#Add missing IDSize, IDReef and IDSpecies columns
#Sometimes the raw db contains the IDSize column, but not Size column
#So it should be adjusted in here
t <-alex %>% left_join(sizes) %>% filter(!is.na(Year)) %>% left_join(reefs) %>% 
  left_join(species)


#Check for missing IDs or Sizes  
missing <- t %>% 
  select(IDReef,Reef) %>% 
  filter(is.na(Reef)) %>% 
  distinct()


#Write the preprocessed db to the respective folder (date) and name accordingly to the Observer column

writexl::write_xlsx(t, "data/raw/2025/OCT-NOV/preprocessed/NAME")


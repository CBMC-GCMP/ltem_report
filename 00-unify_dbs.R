library(tidyverse)
library(readxl)
library(readr)  # for parse_number()

# Path to your folder
path <- "data/raw/2025/OCT-NOV/preprocessed/"

# 1. List all Excel files
files <- list.files(path, pattern = "\\.xlsx$", full.names = TRUE)

# 2. Read ALL columns as text and bind them
ltem_all_raw <- files %>%
  set_names() %>% 
  map_dfr(function(f) {
    read_excel(f, col_types = "text") %>%  # everything as character
      mutate(source_file = basename(f))
  })

# 3. Columns that must be numeric
numeric_cols <- c(
  "Year", "Month", "Day",
  "IDReef", "Depth", "Transect",
  "IDSpecies", "Size", "Quantity", "Area"
)

# 4. Convert the defined numeric columns to numeric
ltem_all <- ltem_all_raw %>%
  mutate(
    across(
      any_of(numeric_cols),
      ~ parse_number(.)  # more robust than as.numeric()
    )
  ) %>% 
  mutate(Observer=ifelse(is.na(Observer), "Alex", Observer)) %>% 
  select(-c(source_file, Comentarios, 20,21))

# Inspect result
glimpse(ltem_all)


write_xlsx(ltem_all, "data/raw/2025/OCT-NOV/ltem_OCT-NOV_2025.xlsx")

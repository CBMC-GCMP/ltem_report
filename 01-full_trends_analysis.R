# Load necessary libraries
library(mgcv)
library(patchwork)
library(tidyverse) # Includes dplyr, tidyr, etc.
library(lme4)
library(ggplot2)
library(scales)
library(purrr) # For map function
library(dplyr)   # Explicitly load for group_by, summarise, etc.
library(tidyr)   # Explicitly load for nest, unnest functions
library(ggridges)
library(dafishr)
library(tidytext)
library(cowplot)

# Create figures directory if it doesn't exist
if (!dir.exists("figures")) {
  dir.create("figures")
}

# Load and prepare data ---------------------------------------------------------------

# Load main datasets
ltem <- readRDS("data/ltem_historic_updated_2025-05-13.RDS") %>% 
  mutate(degree = round(Latitude, 0)) %>% 
  filter(Family != "Carangidae") %>% 
  select(-1)

ltem.pfa <- readRDS("data/ltem_pfa_2021-2024.RDS") %>% 
  mutate(degree = round(Latitude, 0)) %>% 
  select(-1)

# Combine datasets
ltem <- rbind(ltem, ltem.pfa)

# Load sargassum data
sargazo <- readRDS("data/ltem_historic_sargassum_2025-05-12.RDS")

# Source size check function
source("size_check.R")

# Define trophic color palette for consistent use across plots
trophic_colors <- c(
  "Depredadores solitarios" = "#D73027",      # Rojo intenso
  "Depredadores en cardúmenes" = "#FC8D59",   # Naranja
  "Omnívoros en cardúmen" = "#FEE08B",        # Amarillo
  "Herbívoros en cardúmen" = "#91BFDB",       # Azul claro
  "Crípticos solitarios" = "#4575B4",         # Azul medio
  "Planctívoros" = "#313695"                  # Azul oscuro
)

# Fish Analysis ---------------------------------------------------------------

# Prepare fish data
fish_productivity_data <- ltem %>%
  janitor::clean_names() %>%
  filter(label == "PEC") %>%
  mutate(trophic_level = as.numeric(trophic_level)) %>%
  left_join(., comm_sp) %>%
  droplevels() %>%
  left_join(., group_traits) %>%
  mutate(commercial = ifelse(commercial == "yes", "Commercial", "No commercial")) %>%
  filter(region %in% c("Cabo Pulmo", "Corredor",
                      "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  filter(!(region %in% c("Corredor", "Santa Rosalia") & family %in% c("Haemulidae", "Carangidae") & biomass > 3)) %>% 
  filter(!(region == "Santa Rosalia" & year == 2022 & biomass > 2.5)) %>%
  size_check() %>% 
  filter(size_flagged_outlier != "TRUE")

# Calculate biomass at different levels ----------------------------------------

# Function to select consistently monitored sites
select_consistent_sites <- function(data, min_years_threshold = NULL) {
  total_years <- length(unique(data$year))
  message(paste0("Total years in dataset: ", total_years))
  
  reef_monitoring_consistency <- data %>%
    group_by(region, reef) %>%
    summarise(
      years_monitored = n_distinct(year),
      years_coverage = years_monitored / total_years,
      .groups = "drop"
    ) %>%
    arrange(region, desc(years_monitored))
  
  if (is.null(min_years_threshold)) {
    max_coverage_by_region <- reef_monitoring_consistency %>%
      group_by(region) %>%
      summarise(max_coverage = max(years_coverage), .groups = "drop")
    
    min_max_coverage <- min(max_coverage_by_region$max_coverage)
    min_years_count <- ceiling(min_max_coverage * total_years)
    
    message(paste0("Maximum consistent coverage across all regions: ", 
                  round(min_max_coverage * 100, 1), "% (", min_years_count, " years)"))
    min_years_threshold <- min_years_count
  }
  
  consistent_reefs <- reef_monitoring_consistency %>%
    filter(years_monitored >= min_years_threshold)
  
  consistent_sites_count <- consistent_reefs %>%
    group_by(region) %>%
    summarise(site_count = n(), .groups = "drop")
  
  message("Number of consistent sites by region (monitored for at least ", 
         min_years_threshold, " years):")
  print(consistent_sites_count)
  
  filtered_data <- data %>%
    inner_join(consistent_reefs %>% select(region, reef), 
               by = c("region", "reef"))
  
  return(filtered_data)
}

# Calculate biomass at reef level for all regions except Santa Rosalia
biomass_raw_other <- fish_productivity_data %>%
  filter(region != "Santa Rosalia") %>%
  group_by(year, region, protection_level, reef, depth, transect) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  group_by(year, region, protection_level, reef) %>%
  summarise(biomass = mean(biomass), .groups = "drop")

# Calculate biomass at transect level for Santa Rosalia
biomass_raw_santarosalia <- fish_productivity_data %>%
  filter(region == "Santa Rosalia") %>%
  group_by(year, region, protection_level, reef, depth, transect) %>%
  summarise(biomass = sum(biomass), .groups = "drop")

# Combine both datasets
biomass_raw <- bind_rows(
  biomass_raw_other,
  biomass_raw_santarosalia
)

# Apply consistent site selection
balanced_biomass <- select_consistent_sites(biomass_raw)

# Calculate summary statistics by region
biomass_by_region <- balanced_biomass %>%
  group_by(year, region) %>%
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    sd_biomass = sd(biomass, na.rm = TRUE),
    n = n(),
    se_biomass = sd_biomass / sqrt(n),
    .groups = "drop"
  )

# Calculate summary statistics by protection level
biomass_by_protection <- balanced_biomass %>%
  group_by(year, protection_level) %>%
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(protection_level = factor(protection_level, 
                                 levels = c("Prohibited", "Allowed", "Open Area"),
                                 labels = c("Sin Pesca", "Multi-uso", "Sin Protección")))

# Plot biomass trends by region
biomass_by_region %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  ggplot(aes(x = year, y = mean_biomass, color = region)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 4, linetype = 2, col = "firebrick", linewidth = 1) +
  geom_errorbar(
    aes(
      ymin = mean_biomass - se_biomass,
      ymax = mean_biomass + se_biomass
    ),
    width = 0.2
  ) +
  labs(x = "Año", y = "Biomasa Promedio (T/ha)") +
  scale_x_continuous(breaks = seq(1998, 2025, by = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(angle = 0, face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = "Set2")

ggsave("figures/biomass_trends_by_region_2025_2.png", width = 10, height = 6, dpi = 300)

# Plot biomass trends by protection level
ggplot(biomass_by_protection, aes(x = year, y = mean_biomass, color = protection_level)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass, 
                    ymax = mean_biomass + se_biomass),
                width = 0.2) +
  labs(x = "Año", y = "Biomasa Promedio (T/ha)") +
  scale_x_continuous(breaks = seq(1998, 2025, by = 1)) +
  geom_hline(yintercept = 4, linetype = 2, col = "firebrick", linewidth = 1) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(angle = 0, face = "bold"),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = comma) +
  scale_color_brewer(palette = "Set2")

ggsave("figures/biomass_trends_by_protection_2025.png", width = 10, height = 6, dpi = 300)

# Functional Group Analysis ---------------------------------------------------------------

# Calculate biomass by functional group
biomass_f <- fish_productivity_data %>%
  group_by(year, region, protection_level, reef, depth, transect, functional_name) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, functional_name) %>%
  summarise(biomass = mean(biomass), .groups = "drop")

# Apply consistent site selection
balanced_biomass_f <- select_consistent_sites(biomass_f)

# Calculate summary statistics by functional group
biomass_by_functional <- balanced_biomass_f %>%
  group_by(year, region, protection_level, functional_name) %>%
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(protection_level = factor(protection_level, 
                                 levels = c("Prohibited", "Allowed", "Open Area"),
                                 labels = c("Sin Pesca", "Multi-uso", "Sin Protección")))

# Top 10 Species Analysis ---------------------------------------------------------------

# Prepare data for top 10 species analysis
fish.spp <- fish_productivity_data %>%
  group_by(year, region, protection_level, reef, depth, transect, functional_name, species) %>%
  summarise(biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, functional_name, species) %>%
  summarise(biomass = mean(biomass, na.rm = TRUE), .groups = "drop")

# Apply consistent site selection and add period classification
balanced_fish <- select_consistent_sites(fish.spp) %>% 
  mutate(period = case_when(
    year < 2014 ~ "Histórico (1998-2013)",
    year >= 2014 & year <= 2024 ~ "Calentamiento (2014-2024)",
    year == 2025 ~ "2025"
  ))

# Create functional group lookup table
functional.fish <- balanced_fish %>% 
  select(species, functional_name) %>% 
  distinct()

# Three-period version of top 10 species plots
p3 <- balanced_fish %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(degree = ifelse(round(lat, 0) >= 26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  filter(!is.na(functional_name)) %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree != "Debajo de 26°") %>% 
  filter(functional_name != "Pelagic") %>% 
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>%
  mutate(period = factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025"))) %>%
  ungroup() %>% 
  group_by(period) %>%
  mutate(total_biomass = sum(biomass, na.rm = TRUE),
         rel_biomass = (biomass / total_biomass) * 100) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  mutate(functional_name = factor(functional_name, 
                                levels = c("GenPred_solitary", "GenPred_schooling",
                                         "EpiBent_schooling", "Crip_schooling",
                                         "Crip_solitary", "Plank", "Pelagic"),
                                labels = c("Depredadores solitarios",
                                         "Depredadores en cardúmenes",
                                         "Omnívoros en cardúmen",
                                         "Herbívoros en cardúmen",
                                         "Crípticos solitarios",
                                         "Planctívoros",
                                         "Pelágicos"))) %>% 
  filter(functional_name != "Pelágicos") %>% 
  ungroup() %>% 
  mutate(species = reorder_within(species, rel_biomass, period)) %>% 
  ggplot(aes(x = species, y = rel_biomass, fill = functional_name)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5)) +
  coord_flip() +
  facet_wrap(~period, scales = "free_y", nrow = 1) +
  scale_x_reordered() +
  labs(y = "Biomasa Relativa (%)", x = "", fill = "Grupo Funcional", 
       title = "Latitud mayor a 26°") +
  theme_classic() +
  theme(legend.position = "",
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.title.position = "bottom",
        legend.title = element_text(face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = trophic_colors)

# Second plot for latitudes below 26°
p4 <- balanced_fish %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(degree = ifelse(round(lat, 0) >= 26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  filter(!is.na(functional_name)) %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree != "Arriba de 26°") %>% 
  filter(functional_name != "Pelagic") %>% 
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>%
  mutate(period = factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025"))) %>%
  ungroup() %>% 
  group_by(period) %>%
  mutate(total_biomass = sum(biomass, na.rm = TRUE),
         rel_biomass = (biomass / total_biomass) * 100) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  mutate(functional_name = factor(functional_name, 
                                levels = c("GenPred_solitary", "GenPred_schooling",
                                         "EpiBent_schooling", "Crip_schooling",
                                         "Crip_solitary", "Plank", "Pelagic"),
                                labels = c("Depredadores solitarios",
                                         "Depredadores en cardúmenes",
                                         "Omnívoros en cardúmen",
                                         "Herbívoros en cardúmen",
                                         "Crípticos solitarios",
                                         "Planctívoros",
                                         "Pelágicos"))) %>% 
  filter(functional_name != "Pelágicos") %>% 
  ungroup() %>% 
  mutate(species = reorder_within(species, rel_biomass, period)) %>% 
  ggplot(aes(x = species, y = rel_biomass, fill = functional_name)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5)) +
  coord_flip() +
  facet_wrap(~period, scales = "free_y", nrow = 1) +
  scale_x_reordered() +
  labs(y = "Biomasa Relativa (%)", x = "", fill = "Grupo Funcional", 
       title = "Latitud menor a 26°") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.title.position = "bottom",
        legend.title = element_text(face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = trophic_colors)

# Save the three-period version
ggsave("figures/top10_fish_biomass_three_periods.jpg", 
       plot = p3 / p4, 
       width = 15, height = 10, dpi = 300)

# Two-period version ---------------------------------------------------------------

# Prepare data for two-period version
balanced_fish_two_periods <- select_consistent_sites(fish.spp) %>% 
  mutate(period = case_when(
    year < 2025 ~ "Histórico (1998-2024)",
    year == 2025 ~ "2025"
  ))

# Create plots for two-period version
p3_two_periods <- balanced_fish_two_periods %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(degree = ifelse(round(lat, 0) >= 26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  filter(!is.na(functional_name)) %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree != "Debajo de 26°") %>% 
  filter(functional_name != "Pelagic") %>% 
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>%
  mutate(period = factor(period, levels = c("Histórico (1998-2024)", "2025"))) %>%
  ungroup() %>% 
  group_by(period) %>%
  mutate(total_biomass = sum(biomass, na.rm = TRUE),
         rel_biomass = (biomass / total_biomass) * 100) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  mutate(functional_name = factor(functional_name, 
                                levels = c("GenPred_solitary", "GenPred_schooling",
                                         "EpiBent_schooling", "Crip_schooling",
                                         "Crip_solitary", "Plank", "Pelagic"),
                                labels = c("Depredadores solitarios",
                                         "Depredadores en cardúmenes",
                                         "Omnívoros en cardúmen",
                                         "Herbívoros en cardúmen",
                                         "Crípticos solitarios",
                                         "Planctívoros",
                                         "Pelágicos"))) %>% 
  filter(functional_name != "Pelágicos") %>% 
  ungroup() %>% 
  mutate(species = reorder_within(species, rel_biomass, period)) %>% 
  ggplot(aes(x = species, y = rel_biomass, fill = functional_name)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5)) +
  coord_flip() +
  facet_wrap(~period, scales = "free_y", nrow = 1) +
  scale_x_reordered() +
  labs(y = "Biomasa Relativa (%)", x = "", fill = "Grupo Funcional", 
       title = "Latitud mayor a 26°") +
  theme_classic() +
  theme(legend.position = "",
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.title.position = "bottom",
        legend.title = element_text(face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = trophic_colors)

p4_two_periods <- balanced_fish_two_periods %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(degree = ifelse(round(lat, 0) >= 26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  filter(!is.na(functional_name)) %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree != "Arriba de 26°") %>% 
  filter(functional_name != "Pelagic") %>% 
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>%
  mutate(period = factor(period, levels = c("Histórico (1998-2024)", "2025"))) %>%
  ungroup() %>% 
  group_by(period) %>%
  mutate(total_biomass = sum(biomass, na.rm = TRUE),
         rel_biomass = (biomass / total_biomass) * 100) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  mutate(functional_name = factor(functional_name, 
                                levels = c("GenPred_solitary", "GenPred_schooling",
                                         "EpiBent_schooling", "Crip_schooling",
                                         "Crip_solitary", "Plank", "Pelagic"),
                                labels = c("Depredadores solitarios",
                                         "Depredadores en cardúmenes",
                                         "Omnívoros en cardúmen",
                                         "Herbívoros en cardúmen",
                                         "Crípticos solitarios",
                                         "Planctívoros",
                                         "Pelágicos"))) %>% 
  filter(functional_name != "Pelágicos") %>% 
  ungroup() %>% 
  mutate(species = reorder_within(species, rel_biomass, period)) %>% 
  ggplot(aes(x = species, y = rel_biomass, fill = functional_name)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 3.5, 0.5)) +
  coord_flip() +
  facet_wrap(~period, scales = "free_y", nrow = 1) +
  scale_x_reordered() +
  labs(y = "Biomasa Relativa (%)", x = "", fill = "Grupo Funcional", 
       title = "Latitud menor a 26°") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.title.position = "bottom",
        legend.title = element_text(face = "bold", hjust = 0.5)) +
  scale_fill_manual(values = trophic_colors)

# Save the two-period version
ggsave("figures/top10_fish_biomass_two_periods.jpg", 
       plot = p3_two_periods / p4_two_periods, 
       width = 15, height = 10, dpi = 300)

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
# INV ---------------------------------------------------------------------

# Load main invertebrate data -----
invertebrate_data <- ltem %>%
  janitor::clean_names() %>%
  filter(label == "INV") %>%  # Focus on invertebrates instead of fish
  filter(region %in% c("Cabo Pulmo", "Corredor",
                       "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia"))



## Filter data for Holaxonia and Scleractinia ----
coral_data <- invertebrate_data %>%
  filter(taxa3 %in% c("Holaxonia", "Scleractinia")) %>%
  # Handle quantity - assuming 'quantity' is the column name
  # If it's different, adjust accordingly (biomass, density, etc.)
  group_by(year, region, protection_level, reef, depth, transect, taxa3) %>%
  summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, taxa3) %>%
  summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop")

# Apply the function to get the consistently monitored sites
balanced_coral_data <- select_consistent_sites(coral_data)

# Calculate mean quantity by year, region, and order
coral_by_region <- balanced_coral_data %>%
  dplyr::group_by(year, region, taxa3) %>%
  dplyr::summarise(
    mean_quantity = mean(quantity, na.rm = TRUE),
    se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Calculate mean quantity by year, protection level, and order
coral_by_protection <- balanced_coral_data %>%
  dplyr::group_by(year, protection_level, taxa3) %>%
  dplyr::summarise(
    mean_quantity = mean(quantity, na.rm = TRUE),
    se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Plot trends by region

coral_by_region %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(region=factor(region, levels = c("Santa Rosalia",
                                          "San Basilio",
                                          "Loreto",
                                          "Corredor",
                                          "La Paz", 
                                          "La Ventana",
                                          "Cabo Pulmo"))) %>% 
  filter(!region%in% c("Cabo Pulmo", "La Ventana")) %>% 
ggplot(
                            aes(x = year, y = mean_quantity, color = taxa3)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_quantity - se_quantity,
      ymax = mean_quantity + se_quantity
    ),
    width = 0.2
  ) +
  scale_x_continuous(breaks = seq(1998, 2025, by = 1)) +
  labs(x = "Año", y = "Abundancia Promedio",
       # title = "Holaxonia vs Scleractinia quantity by Region"
       ) +
  # scale_x_continuous(breaks = seq(min(coral_by_region$year), max(coral_by_region$year), by = 2)) +
  theme_bw() +
  facet_wrap(~region, ncol=1)+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color=NA, fill=NA),
        strip.text = element_text(face="bold", size=10),
        axis.text.x=element_text(angle=90),
        axis.title = element_text(face="bold")) +

  scale_y_continuous(labels = comma) +
  scale_color_viridis_d(direction = -1)
# theme(legend.position = "bottom")
ggsave("figures/inv_trends_by_region_2025_region.png",  width = 10, height = 12, dpi = 300)

# Plot trends by protection level
coral_protection_plot <- ggplot(coral_by_protection, 
                                aes(x = year, y = mean_quantity, 
                                    color = protection_level, linetype = taxa3)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_quantity - se_quantity, 
                    ymax = mean_quantity + se_quantity),
                width = 0.2) +
  labs(x = "Año", y = "Abundancia Promedio", 
       # title = "Holaxonia vs Scleractinia quantity by Protection Level"
       ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(1998, 2025, by = 1)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.x=element_text(angle=90),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(labels = comma) +
  scale_color_viridis_d()
# theme(legend.position = "bottom")

# Combine plots using patchwork
coral_combined_plot <- coral_region_plot + coral_protection_plot + plot_layout(ncol = 1)
print(coral_combined_plot)

# Fit GAM models to assess trends
# By region and order
coral_by_region$region_factor <- as.factor(coral_by_region$region)
coral_by_region$taxa3_factor <- as.factor(coral_by_region$taxa3)

gam_region_coral <- gam(mean_quantity ~ s(year, by = region_factor:taxa3_factor) + 
                          region_factor * taxa3_factor, 
                        data = coral_by_region, 
                        family = gaussian())
summary(gam_region_coral)

# By protection level and order
coral_by_protection$protection_factor <- as.factor(coral_by_protection$protection_level)
coral_by_protection$taxa3_factor <- as.factor(coral_by_protection$taxa3)

gam_protection_coral <- gam(mean_quantity ~ s(year, by = protection_factor:taxa3_factor) + 
                              protection_factor * taxa3_factor, 
                            data = coral_by_protection, 
                            family = gaussian())
summary(gam_protection_coral)









# Save the plots

ggsave("figures/inv_trends_by_protection_2025.png", coral_protection_plot, width = 10, height = 6, dpi = 300)





## Asteroidea and Echinoidea ----------

# Filter data for Asteroidea and Echinoidea
echinoderm_data <- invertebrate_data %>%
  filter(taxa2 %in% c("Asteroidea", "Echinoidea")) %>%
  # Calculate richness at the reef level
  group_by(year, region, protection_level, reef, depth, transect, taxa2) %>%

  summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, taxa2) %>%
  summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop")
# Apply the function to get the consistently monitored sites
balanced_echinoderm_data <- select_consistent_sites(echinoderm_data)

# Calculate mean richness by year, region, and class
echinoderm_by_region <- balanced_echinoderm_data %>%
  dplyr::group_by(year, region, taxa2) %>%
  dplyr::summarise(
    mean_quantity = mean(quantity, na.rm = TRUE),
    se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Calculate mean quantity by year, protection level, and class
echinoderm_by_protection <- balanced_echinoderm_data %>%
  dplyr::group_by(year, protection_level, taxa2) %>%
  dplyr::summarise(
    mean_quantity = mean(quantity, na.rm = TRUE),
    se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Plot trends by region

echinoderm_by_region %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(region=factor(region, levels = c("Santa Rosalia",
                                          "San Basilio",
                                          "Loreto",
                                          "Corredor",
                                          "La Paz", 
                                          "La Ventana",
                                          "Cabo Pulmo"))) %>% 
  filter(!region%in% c("Cabo Pulmo", "La Ventana")) %>% 
  mutate(log_quan=log1p(mean_quantity),
         se_log=log1p(se_quantity)) %>% 
ggplot( 
                                 aes(x = year, y = log_quan, color = taxa2)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = log_quan - se_log,
      ymax = log_quan + se_log
    ),
    width = 0.2
  ) +
  labs(x = "Año", y = "Abundancia Promedio (Escala logarítmica)") +
  scale_x_continuous(breaks = seq(min(echinoderm_by_region$year), 
                                  max(echinoderm_by_region$year), by = 1)) +
  theme_bw() +
   facet_wrap(~region, ncol=1)+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(color=NA, fill=NA),
        strip.text = element_text(face="bold", size=10),
        axis.text.x=element_text(angle=90),
        axis.title = element_text(face="bold")) +
  scale_y_continuous(labels = comma) +
  scale_color_viridis_d()

 ggsave("figures/inv_trends_ech_ast_abund_2025_log.png",  width = 10, height = 12, dpi = 300)


 

# INv bar plots -----
 

# Filter data for Asteroidea and Echinoidea
inv.abund  <- invertebrate_data %>%
   mutate(taxa=case_when(taxa3== "Holaxonia"~"Holaxonia",
                         taxa3=="Scleractinia"~ "Scleractinia",
                         taxa2=="Asteroidea"~"Asteroidea",
                         taxa2=="Echinoidea"~"Echinoidea")) %>% 
   filter(taxa %in% c("Asteroidea", "Echinoidea", "Scleractinia", "Holaxonia")) %>%
   # Calculate richness at the reef level
   group_by(year, region, protection_level, reef, depth, transect, taxa) %>%
   
   summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
   group_by(year, region, protection_level, reef, taxa) %>%
   summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop")
# Apply the function to get the consistently monitored sites
 balanced_inv <- select_consistent_sites(inv.abund)
 
# 
invby_region <- balanced_inv %>%
   dplyr::group_by(year, region, taxa) %>%
   dplyr::summarise(
     mean_quantity = mean(quantity, na.rm = TRUE),
     se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
     .groups = "drop"
   )


absolute.bars.inv <- invby_region %>% 
  filter(year==2025) %>% 
  # filter(region!="La Ventana", "Cabo Pulmo" ) %>% 
  
# mutate(category=ifelse(year=2025, "2025", "Histórico"))

  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=mean_quantity, fill=taxa, col=taxa))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
  labs(y="Abundancia Promedio", x="Región", fill="Grupo Taxonómico")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "",
        legend.title = element_text(hjust=0.5, face="bold"))


relative.inv <- invby_region %>% 
  # filter(region!="La Ventana") %>% 
  filter(year==2025) %>% 


  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=mean_quantity, fill=taxa, col=taxa))+
  geom_bar(stat="identity", position="fill")+
  scale_y_continuous(labels = label_percent(scale = 100)) + 
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
  labs(y="Abundancia Promedio (%)", x="Región", fill="Grupo Taxonómico")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        legend.position = "left",
        axis.title = element_text(face="bold", size=10),
        legend.title = element_text(hjust=0.5, face="bold"))
(biomass.bars <- absolute.bars.inv + relative.inv + plot_layout(ncol = 2))

ggsave("figures/invabund_bars_2025_region.png", biomass.bars, width = 12, height = 6, dpi = 300)


### Top 10 Inv Abundances ------

# Filter data for Asteroidea and Echinoidea
inv.spp <- invertebrate_data %>%
  mutate(taxa=case_when(taxa3== "Holaxonia"~"Holaxonia",
                        taxa3=="Scleractinia"~ "Scleractinia",
                        taxa2=="Asteroidea"~"Asteroidea",
                        taxa2=="Echinoidea"~"Echinoidea")) %>% 
  filter(taxa %in% c("Asteroidea", "Echinoidea", "Scleractinia", "Holaxonia")) %>%
  mutate(species=case_when(str_detect(species, "Muricea")~"Muricea sp",
                           T~species)) %>% 
  # Calculate richness at the reef level
  group_by(year, region, protection_level, reef, depth, transect, taxa, species) %>%
  
  summarise(quantity = sum(quantity, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, taxa, species) %>%
  summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop")
# Apply the function to get the consistently monitored sites
balanced_inv <- select_consistent_sites(inv.spp) %>% 
  mutate(period=case_when( year< 2014 ~"Histórico (1998-2013)",
                         year>= 2014 & year<=2024 ~"Calentamiento (2014-2024)",
                         year ==2025 ~"2025"))








taxa.inv <- taxa.inv <- balanced_inv %>% 
  select(taxa, species) %>% 
  distinct()


 p3 <- balanced_inv %>% 
  # filter(year>2021) %>% 
   left_join(lats %>% janitor::clean_names()) %>% 
   mutate(degree=ifelse(round(lat,0)>=26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree!="Arriba de 26°") %>% 
   mutate(period=factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025")),
          species=ifelse(str_detect(species, "Pocillopora"), "Pocillopora sp", species)
          ) %>%

   
  group_by(period, taxa, species) %>% 
  summarise(quantity = mean(quantity, na.rm=T)) %>%
  # filter(MPA=="MPA") %>%
  # mutate(year=factor(year)) %>% 
  group_by(period) %>%
  top_n(10, quantity) %>%
  ungroup() %>% 
  mutate(species= reorder_within(species, quantity, period)) %>% 
  ggplot(aes(x=species, y = quantity, fill=taxa)) +
  geom_col()+
  coord_flip()+
  facet_wrap(~period, scales="free_y", nrow=1)+
  scale_x_reordered() +
  # scale_color_material_d() +
  # scale_fill_manual(values = c("firebrick", "darkgreen"))+
  labs(y = "Abundancia Promedio", x="", fill="Grupo Taxonómico",  title= "Latitud menor a 26°") +
  theme_classic() +
   theme(legend.position = "") +
   guides(colour = "none")+
   theme(axis.text.y = element_text(face= "italic"),
         axis.title = element_text(face="bold"),
         plot.title = element_text(hjust=0.5, face="bold", size=15),
         strip.background = element_blank(),  # Removes the background
         strip.text = element_text(face = "bold", hjust = 0.5, size=12),
         legend.title.position = "bottom",
         legend.title = element_text(face="bold", hjust=0.5))+
   scale_fill_viridis_d()


 
 

 p4 <- balanced_inv %>% 
   # filter(year>2021) %>% 
   left_join(lats %>% janitor::clean_names()) %>% 
   mutate(degree=ifelse(round(lat,0)>=26, "Arriba de 26°", "Debajo de 26°")) %>%  
   filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree!="Debajo de 26°") %>% 
   mutate(period=factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025")),
          species=ifelse(str_detect(species, "Pocillopora"), "Pocillopora sp", species)
   ) %>%
   
   
   group_by(period, taxa, species) %>% 
   summarise(quantity = mean(quantity, na.rm=T)) %>%
   # filter(MPA=="MPA") %>%
   # mutate(year=factor(year)) %>% 
   group_by(period) %>%
   top_n(10, quantity) %>%
   ungroup() %>% 
   mutate(species= reorder_within(species, quantity, period)) %>% 
   ggplot(aes(x=species, y = quantity, fill=taxa)) +
   geom_col()+
   coord_flip()+
   facet_wrap(~period, scales="free_y", nrow=1)+
   scale_x_reordered() +
   # scale_color_material_d() +
   # scale_fill_manual(values = c("firebrick", "darkgreen"))+
   labs(y = "", x="", fill="Grupo Taxonómico", title= "Latitud mayor a 26°") +
   theme_classic() +
   theme(legend.position = "bottom") +
   guides(colour = "none")+
   theme(axis.text.y = element_text(face= "italic"),
         axis.title = element_text(face="bold"),
         plot.title = element_text(hjust=0.5, face="bold", size=15),
         strip.background = element_blank(),  # Removes the background
         strip.text = element_text(face = "bold", hjust = 0.5, size=12),
         legend.title.position = "bottom",
         
         legend.title = element_text(face="bold", hjust=0.5))+
   scale_fill_viridis_d()


plot_grid(
  p4 + theme(legend.position = "none"),  # Remove legend from p1
  p3 + theme(legend.position = "bottom"),  # Remove legend from p2
  labels = c("A", "B"),   
  
  # rel_heights = c(3, 0.1),
  # Add labels A and B
  ncol = 1                               # Arrange plots in one column
)
ggsave("figures/top10_inv_abund_hist.jpg", width = 12, height = 12, dpi = 300)

# Diversity Indexes -------------------------------------------------------
richness <-  ltem %>%
   janitor::clean_names() %>% 
   filter(taxa2 %in% c("Asteroidea", "Echinoidea")|taxa3 %in% c("Holaxonia", "Scleractinia")|label=="PEC" ) %>%

 group_by(year, region, protection_level, reef, depth, transect, label) %>%
   summarise(richness = n_distinct(species), .groups = "drop") %>%
   group_by(year, region, protection_level, reef, label) %>%
   summarise(richness = mean(richness), .groups = "drop")
 
 # Apply the function to get the consistently monitored sites
 balanced_richness <- select_consistent_sites(richness)
 
 
 
 # Calculate mean richness by year, region, and class
 richness_by_region <- balanced_richness %>%
   dplyr::group_by(year, region, label) %>%
   dplyr::summarise(
     mean_richness = mean(richness, na.rm = TRUE),
     se_richness = sd(richness, na.rm = TRUE) / sqrt(dplyr::n()),
     .groups = "drop"
   )
 
 
richness_by_region %>% 
  filter(region %in% c("Corredor",
                       "La Paz",  "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  mutate(category=ifelse(year<2025, "Histórico", "2025"),
         category=factor(category, levels=c("Histórico", "2025"))) %>% 
  group_by(category, region, label) %>% 
  summarise(richness=mean(mean_richness, na.rm=T)) %>% 
  filter(label=="INV") %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, richness), y=richness, fill= region)) +
  geom_bar(stat="identity")+
  facet_wrap(~category, ncol=2)+
  labs(x="Región", y="Riqueza de especies (Invertebrados)")+
  coord_flip()+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(axis.title = element_text(face="bold", size=10),
        legend.position = "")



## Shannon Index -----
library(vegan)
library(dplyr)
library(tidyr)
library(vegan)
library(tibble)

shannon_transect_hist <- ltem %>%
  filter(!is.na(Quantity), !is.na(Species)) %>%
  group_by(Region, Reef, Transect, Species) %>%
  summarise(total_quantity = sum(Quantity, na.rm = TRUE), .groups = "drop") %>%
  unite(transect_id, Region, Reef, Transect, remove = FALSE) %>%
  pivot_wider(names_from = Species, values_from = total_quantity, values_fill = 0) %>%
  # Asegura que solo las columnas numéricas vayan a vegan::diversity
  {
    data_wide <- .
    species_matrix <- data_wide %>%
      select(where(is.numeric)) %>%
      as.data.frame()
    
    shannon_values <- vegan::diversity(species_matrix, index = "shannon")
    
    tibble::tibble(transect_id = data_wide$transect_id,
                   Region = data_wide$Region,
                   Reef = data_wide$Reef,
                   Transect = data_wide$Transect,
                   Shannon_Index = shannon_values)
  }

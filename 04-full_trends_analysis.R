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


#Custom size check function
size_check <- function(data, size_col = "size", species_col = "species") {
  data %>%
    group_by(across(all_of(species_col))) %>%
    mutate(
      q1 = quantile(.data[[size_col]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[size_col]], 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lower_bound = q1 - 1.5 * iqr,
      upper_bound = q3 + 1.5 * iqr,
      size_flagged_outlier = .data[[size_col]] < lower_bound | .data[[size_col]] > upper_bound
    ) %>%
    ungroup()
}



# Load and prepare data ---------------------------------------------------------------

#Load comercial spp
comm_sp <- read.csv("data/commercial_species.csv")

group_traits <- read.csv("data/cluster_to_create_traits.csv") %>% 
  janitor::clean_names()
# Load main datasets
ltem <- readRDS("outputs/historical/ltem_historic_updated_2025-11-20.RDS") %>% 
  mutate(degree = round(Latitude, 0))


ltem.pfa <- readRDS("data/ltem_pfa_2021-2024.RDS") %>% 
  mutate(degree = round(Latitude, 0)) %>% 
  select(-1)

# Combine datasets
ltem <- rbind(ltem, ltem.pfa)


lats <- ltem %>% 
  filter(Region %in% c( "Cabo Pulmo", "Corredor",
                        "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  select(Region, Latitude) %>% 
  group_by(Region) %>% 
  summarise(lat=mean(Latitude, na.rm=T))

# Load sargassum data
# sargazo <- readRDS("data/ltem_historic_sargassum_2025-05-12.RDS")



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
  filter(Family != "Carangidae") %>% 
  janitor::clean_names() %>%
  filter(label == "PEC") %>%
  mutate(trophic_level = as.numeric(trophic_level)) %>%
  left_join(., comm_sp) %>%
  droplevels() %>%
  left_join(., group_traits) %>%
  mutate(commercial = ifelse(commercial == "yes", "Commercial", "No commercial")) %>%
  filter(region %in% c("Cabo Pulmo", "Corredor",
                      "La Paz",  "Loreto" )) %>% 
  filter(! (region=="Corredor" & family %in% c("Haemulidae", "Carangidae") & biomass > 3)) %>% 
  size_check() %>%
  filter(size_flagged_outlier!="TRUE")

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


# Calculate biomass at the reef level
biomass_raw <- fish_productivity_data %>%
  group_by(year, region, protection_level, reef, depth, transect) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  group_by(year, region, protection_level, reef) %>%
  summarise(biomass = mean(biomass), .groups = "drop")



# Apply consistent site selection
balanced_biomass <- select_consistent_sites(biomass_raw, 5)

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
  # filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  ggplot(aes(x = year, y = mean_biomass, color = region)) +
  geom_point(size = 2) +
  geom_line() +
  # geom_hline(yintercept = 4, linetype = 2, col = "firebrick", linewidth = 1) +
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
  scale_color_brewer(palette = "Set1")

ggsave("figures/figure1_biomass_trend_per_region.png", width = 10, height = 6, dpi = 300)



# Functional Group Analysis ---------------------------------------------------------------

# Calculate biomass by functional group
biomass_f <- fish_productivity_data %>%
  group_by(year, region, protection_level, reef, depth, transect, functional_name) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, functional_name) %>%
  summarise(biomass = mean(biomass), .groups = "drop")

# Apply consistent site selection
balanced_biomass_f <- select_consistent_sites(biomass_f, 5)

# Calculate summary statistics by functional group
biomass_by_functional <- balanced_biomass_f %>%
  group_by(year, region, protection_level, reef, functional_name) %>%
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(protection_level = factor(protection_level, 
                                 levels = c("Prohibited", "Allowed", "Open Area"),
                                 labels = c("Sin Pesca", "Multi-uso", "Sin Protección"))) %>% 
  mutate( functional_name=factor(functional_name, levels=c("GenPred_solitary",
                                                           "GenPred_schooling",
                                                           "EpiBent_schooling",
                                                           "Crip_schooling",
                                                           "Crip_solitary",
                                                           "Plank",
                                                           "Pelagic"),
                                 labels=c("Depredadores solitarios",
                                          "Depredadores en cardúmenes",
                                          "Omnívoros en cardúmen",
                                          "Herbívoros en cardúmen",
                                          "Crípticos solitarios",
                                          
                                          "Planctívoros",
                                          "Pelágicos"))) %>% 
  filter(functional_name!="Pelágicos", !is.na(functional_name))



# Define trophic color palette for consistent use across plots
trophic_colors <- c(
  "Depredadores solitarios" = "#D73027",      # Rojo intenso
  "Depredadores en cardúmenes" = "#FC8D59",   # Naranja
  "Omnívoros en cardúmen" = "#FEE08B",        # Amarillo
  "Herbívoros en cardúmen" = "#91BFDB",       # Azul claro
  "Crípticos solitarios" = "#4575B4",         # Azul medio
  "Planctívoros" = "#313695"                  # Azul oscuro
)



## Functional groups time series------------

### Full region comparison ----
biomass_by_functional %>% 
  group_by(year, region, functional_name) %>% 
  summarise(biomass=mean(mean_biomass, na.rm=T)) %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(region = reorder(region, -lat)) %>% 
  ggplot(aes(x=year, y=biomass, fill=functional_name, col=functional_name))+
  geom_bar(stat="identity")+
  facet_wrap(~region,ncol=1)+
  scale_fill_manual(values=trophic_colors)+ scale_color_manual(values = trophic_colors, guide = "none")+
  theme_bw()+
  labs(x="Año", y="Biomasa (ton/ha)", fill="Grupo Funcional")+
  scale_x_continuous(breaks = seq(min(biomass_by_functional$year  ), max(biomass_by_functional$year), 1))+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        
        legend.position = "bottom",
        plot.title = element_text(face="bold", size=14, hjust=0.5),
        legend.title = element_text(hjust=0.5, face="bold"), 
        legend.title.position = "bottom")+
  ylim(0, 4)


### Region specific ----
plot_biomass_region <- function(region_name) {
  
  biomass_by_functional %>% 
    filter(region == region_name) %>% 
    group_by(year, region, functional_name) %>% 
    summarise(biomass = mean(mean_biomass, na.rm = TRUE), .groups = "drop") %>% 
    left_join(lats %>% janitor::clean_names(), by = "region") %>% 
    ggplot(aes(x = year, 
               y = biomass, 
               fill = functional_name, 
               col = functional_name)) +
    
    geom_bar(stat = "identity") +
    
    scale_fill_manual(values = trophic_colors) +
    scale_color_manual(values = trophic_colors, guide = "none") +
    
    theme_bw() +
    labs(
      title = region_name,
      x = "Año",
      y = "Biomasa (ton/ha)",
      fill = "Grupo Funcional"
    ) +
    
    scale_x_continuous(
      breaks = seq(min(biomass_by_functional$year),
                   max(biomass_by_functional$year),
                   1)
    ) +
    
    theme(
      axis.text.x = element_text(angle = 0),
      axis.title = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.title = element_text(hjust = 0.5, face = "bold"),
      legend.title.position = "bottom"
    ) +
    
    ylim(0, 4)
}

plot_biomass_region("Cabo Pulmo")

ggsave("figures/figure2_biomass_bars_per_region_CaboPulmo.png", width = 10, height = 6, dpi = 300)




## Functional groups reef comparison historical vs current ----
plot_biomass_bars <- function(region_name) {
  
  # --- helper function to summarize & plot ---
  make_plot <- function(data, title_text, show_legend = FALSE) {
    
    data %>%
      group_by(functional_name) %>%
      summarise(
        biomass = mean(mean_biomass, na.rm = TRUE),
        sd      = sd(mean_biomass,  na.rm = TRUE),
        n       = sum(!is.na(mean_biomass)),
        se      = sd / sqrt(n),
        .groups = "drop"
      ) %>%
      mutate(functional_name = fct_rev(functional_name)) %>%
      ggplot(aes(x = functional_name,
                 y = biomass,
                 fill = functional_name)) +
      
      geom_col(color = NA) +
      geom_errorbar(aes(ymin = biomass - se, ymax = biomass + se),
                    width = 0.2,
                    linewidth = 0.6,
                    color = "black") +
      
      scale_fill_manual(values = trophic_colors) +
      
      labs(
        title = title_text,
        x = "Grupo Funcional",
        y = "Biomasa (ton/ha)",
        fill = "Grupo Funcional"
      ) +
      
      coord_flip() +
      theme_bw() +
      theme(
        axis.title = element_text(face = "bold", size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = ifelse(show_legend, "right", "none"),
        legend.title = element_text(face = "bold")
      ) +
      ylim(0, 1)
  }
  
  # --- CURRENT YEAR PLOT ---
  current_data <- biomass_by_functional %>%
    filter(year == max(ltem$Year), region == region_name)
  
  p_current <- make_plot(current_data, title_text = max(ltem$Year), show_legend = FALSE)
  
  # --- HISTORICAL PLOT ---
  historical_data <- biomass_by_functional %>%
    filter(region == region_name)
  
  p_hist <- make_plot(historical_data, title_text = "Histórico", show_legend = TRUE)
  
  # --- combined output ---
  p_hist + p_current + plot_layout(ncol = 2)
}


plot_biomass_bars("Cabo Pulmo")


ggsave("figures/figure3_biomass_bars_functiona_CaboPulmo.png", width = 10, height = 6, dpi = 300)




## Full regions analysis functional groups historic vs current ----------


### Current year -----
abs.biom.2025 <- biomass_by_functional %>% 
  filter(year==max(ltem$Year)) %>% 
  # filter(region!="La Ventana", ) %>% 
  

  group_by(region, functional_name) %>% 
  summarise(biomass=mean(mean_biomass, na.rm=T)) %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  
  ggplot(aes(x=reorder(region, -lat), y=biomass, fill=functional_name, col=functional_name))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=trophic_colors)+ scale_color_manual(values = trophic_colors, guide = "none")+
  labs(y="", x="Región", fill="Grupo Funcional", title="2025")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "",
        plot.title = element_text(face="bold", size=14, hjust=0.5),
        legend.title = element_text(hjust=0.5, face="bold"))+
  ylim(0,2.25)



### Historical period ----

hist.biom.abs <- biomass_by_functional %>% 
  # filter(year==2025) %>% 
  # filter(!region %in% c("La Ventana","Cabo Pulmo" ) )%>% 

  group_by(region, functional_name) %>% 
  summarise(biomass=mean(mean_biomass, na.rm=T)) %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=biomass, fill=functional_name, col=functional_name))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=trophic_colors)+ scale_color_manual(values = trophic_colors, guide = "none")+
  labs(y="Biomasa Promedio (T/ha)", x="Región", fill="Grupo Funcional", title="Histórico")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "right",
        plot.title = element_text(face="bold", size=14, hjust=0.5),
        legend.title = element_text(hjust=0.5, face="bold"))+
  ylim(0,2.25)


(biomass.bars <- hist.biom.abs + abs.biom.2025 + plot_layout(ncol = 2))





ggsave("figures/figure3_biomass_bars_hist-vs-current.png", biomass.bars, width = 12, height = 6, dpi = 300)


























# NRSI -----
nrsi_bal <- fish_productivity_data %>% select_consistent_sites(5)


consistent_reefs <- nrsi_bal %>% 
  filter(year==max(nrsi_bal$year)) %>% 
  select(region,reef) %>% 
  distinct()


nrsi  <- fish_productivity_data %>% 
  filter(reef %in% consistent_reefs$reef) %>% 
  
  mutate(index_levels = case_when(
    str_detect(trophic_level_f, "4-4.5") ~ "UTL", 
    str_detect(trophic_level_f, "2-2.5") ~ "LTL", 
    TRUE ~ "CTL"
  ), 
  # reef=ifelse(reef=="CORONADO_PUNTA_LOBOS", "CORONADO_LAJAS", reef)
  # protection = factor(protection, levels = c("Not Protected", "Lightly Protected", "Fish Refuge", "Fully Protected"))
  ) %>%
  mutate(degree = round(latitude)) %>% 
  group_by(year,degree,  region,reef,  latitude, longitude, depth2, index_levels, transect) %>% 
  summarise(biomass = sum(biomass, na.rm = TRUE)) %>% 
  group_by(year,degree, region, reef,  index_levels,  latitude, longitude) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>% 
  group_by(year,   degree,region, reef) %>% 
  mutate(rel_biomass = (biomass / sum(biomass, na.rm = TRUE)) * 100) %>% 
  dplyr::select(-biomass) %>%
  pivot_wider(names_from = index_levels, values_from = rel_biomass) %>%
  mutate(
    UTL = if_else(is.na(UTL), 0, UTL),  # Replace NA with 0 for UTL calculations
    LTL = if_else(is.na(LTL), 0, LTL),  # Replace NA with 0 for LTL calculations
    CTL = if_else(is.na(CTL), 0, CTL),  # Replace NA with 0 for CTL calculations
    nrsi = case_when(
      LTL > UTL + CTL ~ UTL / (UTL + CTL),  # Condition 1: Use only UTL if LTL is greater than UTL + CTL
      TRUE ~ (UTL + LTL - CTL) / (UTL + LTL + CTL)  # Condition 2: Compute as normal otherwise
    )
  ) 



resample_mean <- function(df, n = 100) {
  replicate(n, {
    sampled <- df[sample(nrow(df), replace = TRUE), ]
    mean(sampled$nrsi, na.rm = TRUE)
  })
}


plot_nrsi_region <- function(region_name, n_boot = 100) {
  
  # ---- LATITUDES ----
  lats <- ltem %>% 
    filter(Region %in% c("Cabo Pulmo","Corredor","La Paz","La Ventana",
                         "Loreto","San Basilio","Santa Rosalia")) %>% 
    group_by(Region, Reef) %>% 
    summarise(lat = mean(Latitude, na.rm = TRUE), .groups = "drop") %>% 
    janitor::clean_names()
  
  # ---- HISTORICAL BOOTSTRAPPING ----
  resampled_hist <- nrsi %>%
    group_by(reef) %>%
    nest() %>%
    mutate(
      ResampledMeans = map(data, ~ resample_mean(.x, n = n_boot)),
      Median = map_dbl(ResampledMeans, mean)
    ) %>%
    dplyr::select(reef, ResampledMeans, Median)
  
  hist_dat <- resampled_hist %>%
    unnest(ResampledMeans) %>% 
    left_join(lats, by = c("reef"))
  
  # ---- CURRENT YEAR BOOTSTRAPPING ----
  max_year <- max(nrsi_bal$year, na.rm = TRUE)
  
  resampled_current <- nrsi %>%
    dplyr::filter(year == max_year) %>%
    group_by(reef) %>%
    nest() %>%
    mutate(
      ResampledMeans = map(data, ~ resample_mean(.x, n = n_boot)),
      Median = map_dbl(ResampledMeans, mean)
    ) %>%
    dplyr::select(reef, ResampledMeans, Median)
  
  current_dat <- resampled_current %>%
    unnest(ResampledMeans) %>% 
    left_join(lats, by = c("reef"))
  
  # ---- FILTER REGION ----
  hist_region <- hist_dat %>% dplyr::filter(region == region_name)
  current_region <- current_dat %>% dplyr::filter(region == region_name)
  
  # ---- PLOTS (igual que antes) ----
  p_hist <- ggplot(hist_region,
                   aes(y = reorder(reef, lat),
                       x = ResampledMeans,
                       fill = stat(x))) +
    geom_density_ridges_gradient(scale = 1) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "", y = "", title = "Histórico") +
    scale_fill_gradientn(
      colors = c("firebrick", "orange", "gray", "grey", "green"),
      values = scales::rescale(c(-1, -0.75, 0, 0.5, 0.75, 1))
    ) +
    xlim(-1, 1) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    )
  
  p_current <- ggplot(current_region,
                      aes(y = reorder(reef, lat),
                          x = ResampledMeans,
                          fill = stat(x))) +
    geom_density_ridges_gradient(scale = 1) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "NRSI", y = "", title = max_year) +
    scale_fill_gradientn(
      colors = c("firebrick", "orange", "gray", "grey", "grey"),
      values = scales::rescale(c(-1, -0.75, 0, 0.5, 0.75, 1))
    ) +
    xlim(-1, 1) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  p_hist + p_current + patchwork::plot_layout(ncol = 1)
}



plot_nrsi_region("Cabo Pulmo")
ggsave("figures/figure4_NRSI_region_CaboPulmo.png", width = 10, height = 6, dpi = 300)



# INV ---------------------------------------------------------------------

# Filter data for Asteroidea and Echinoidea
inv.abund  <- invertebrate_data %>%
  mutate(taxa=case_when(taxa3== "Holaxonia"~"Holaxonia",
                        taxa3=="Scleractinia"~ "Scleractinia",
                        taxa2=="Asteroidea"~"Asteroidea",
                        taxa2=="Echinoidea"~"Echinoidea")) %>% 
  filter(taxa %in% c("Asteroidea", "Echinoidea", "Scleractinia", "Holaxonia")) %>%
  # Calculate richness at the reef level
  group_by(year, region, protection_level, reef, depth, transect, taxa) %>%
  
  summarise(quantity = sum(quantity, na.rm = TRUE),
            richness=n_distinct(species),
            .groups = "drop") %>%
  group_by(year, region, protection_level, reef, taxa) %>%
  summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop",
            richness=mean(richness))


# Apply the function to get the consistently monitored sites
balanced_inv <- select_consistent_sites(inv.abund, 5)



## Abundance vs Richness -----


# Abundance by region
abundby_region <- balanced_inv %>%
  filter(year==max(balanced_inv$year)) %>% 
  dplyr::group_by( region, taxa) %>%
  dplyr::summarise(
    mean_quantity = mean(quantity, na.rm = TRUE),
    se_quantity = sd(quantity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )


# Richness by region
richnessby_region <- balanced_inv %>%
  filter(year==max(balanced_inv$year)) %>% 
  dplyr::group_by( region, taxa) %>%

  dplyr::summarise(
   richness =mean(richness, na.rm = TRUE),
    se_quantity = sd(richness, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )





lats <- ltem %>% 
  filter(Region %in% c( "Cabo Pulmo", "Corredor",
                        "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  select(Region, Latitude) %>% 
  group_by(Region) %>% 
  summarise(lat=mean(Latitude, na.rm=T))








abundance <- abundby_region %>% 
  # filter(year==max(balanced_inv$year)) %>% 
  
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=mean_quantity, fill=taxa, col=taxa))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
  labs(y="Abundancia Promedio", x="Región", fill="Grupo Taxonómico")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "right",
        legend.title = element_text(hjust=0.5, face="bold"))


richness <- richnessby_region %>% 
  # filter(year==max(balanced_inv$year)) %>% 
  
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=richness, fill=taxa, col=taxa))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
  labs(y="Riqueza de especies", x="Región", fill="Grupo Taxonómico")+
  theme_bw()+ 
  
  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "",
        legend.title = element_text(hjust=0.5, face="bold"))

(inv.bars <- abundance + richness + plot_layout(ncol = 2))


ggsave("figures/figure4_inv_abund-rich_bars_2025_region.png", inv.bars, width = 12, height = 6, dpi = 300)



## Echin vs Asteroidea bar trends -------

plot_inverts_region <- function(region_name) {
  
  balanced_inv %>% 
    filter(region == region_name) %>% 
    filter(taxa %in% c("Asteroidea", "Echinoidea")) %>% 

    group_by(year, taxa) %>% 
    summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop") %>% 
    ggplot(aes(x = year, y = quantity, fill = taxa, color = taxa)) +
    
    geom_bar(stat = "identity") +
    
    scale_x_continuous(
      breaks = seq(min(balanced_inv$year), max(balanced_inv$year), 1)
    ) +
    
    scale_fill_viridis_d() +
    scale_color_viridis_d(guide = "none") +
    
    ylim(0, 100) +
    labs(
      x = "Año",
      y = "Abundancia Promedio",
      fill = ""
    ) +
    
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 90),
      axis.title   = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5, face = "bold")
    )
}


plot_inverts_region("Cabo Pulmo")

ggsave("figures/figure5_echi-ast_bars_CaboPulmo.png",  width = 12, height = 6, dpi = 300)


## Scleractinia vs Holaxonia -------



plot_inverts_region <- function(region_name) {
  
  balanced_inv %>% 
    filter(region == region_name) %>% 
    filter(taxa %in% c("Holaxonia", "Scleractinia")) %>% 
    mutate(taxa=factor(taxa, levels=c("Scleractinia", "Holaxonia"),
                       labels=c("Corales cálidos", "Corales fríos"))) %>% 
    group_by(year, taxa) %>% 
    summarise(quantity = mean(quantity, na.rm = TRUE), .groups = "drop") %>% 
    ggplot(aes(x = year, y = quantity, fill = taxa, color = taxa)) +
    
    geom_bar(stat = "identity") +
    
    scale_x_continuous(
      breaks = seq(min(balanced_inv$year), max(balanced_inv$year), 1)
    ) +
    
    scale_fill_viridis_d(direction=-1) +
    scale_color_viridis_d(direction=-1,guide = "none") +
    
    ylim(0, 150) +
    labs(
      x = "Año",
      y = "Abundancia Promedio",
      fill = ""
    ) +
    
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 90),
      axis.title   = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5, face = "bold")
    )
}


plot_inverts_region("Corredor")

ggsave("figures/figure6_holax-scle_bars_Corredor.png",  width = 12, height = 6, dpi = 300)





  


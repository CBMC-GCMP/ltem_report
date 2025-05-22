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



# Load data ---------------------------------------------------------------


ltem <- readRDS("04-report/outputs/ltem_historic_updated_2025-05-13.RDS") %>% 
  mutate(degree=round(Latitude, 0)) %>% 
  filter(Family!="Carangidae") %>% 
  select(-1)


ltem.pfa <- readRDS("04-report/outputs/ltem_pfa_2021-2024.RDS") %>% 
  mutate(degree=round(Latitude, 0)) %>% 
  select(-1)


ltem <- rbind(ltem, ltem.pfa)

sargazo <- readRDS("04-report/outputs/ltem_historic_sargassum_2025-05-12.RDS")

# Fish --------------------------------------------------------------------


lats <- ltem %>% 
  filter(Region %in% c( "Cabo Pulmo", "Corredor",
                        "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  select(Region, Latitude) %>% 
  group_by(Region) %>% 
  summarise(lat=mean(Latitude, na.rm=T))
## data wrangling ----
comm_sp <- read_csv("data/commercial_species.csv")
group_traits <- read_csv("data/cluster_to_create_traits.csv") %>%
    janitor::clean_names()


fish_productivity_data <-ltem  %>%
  janitor::clean_names() %>%
  filter(label == "PEC") %>%
  mutate(trophic_level = as.numeric(trophic_level)) %>%
  left_join(., comm_sp) |>
  droplevels() %>%
  left_join(., group_traits) %>%
  mutate(commercial = ifelse(commercial == "yes", "Commercial", "No commercial")) %>%
  filter(region %in% c("Cabo Pulmo", "Corredor",
                       "La Paz", "La Ventana", "Loreto", "San Basilio", "Santa Rosalia")) %>% 
  filter(!(region %in% c("Corredor", "Santa Rosalia") & family %in% c("Haemulidae", "Carangidae") & biomass > 3)) %>% 
  filter(!(region=="Santa Rosalia" & year==2022 & biomass > 2.5 ) ) %>%
  size_check() %>% 
  filter(size_flagged_outlier!="TRUE")



## Biomass at reef level ----

# Calculate biomass at the reef level
biomass_raw <- fish_productivity_data %>%
    group_by(year, region, protection_level, reef, depth, transect) %>%
    summarise(biomass = sum(biomass), .groups = "drop") %>%
    group_by(year, region, protection_level, reef) %>%
    summarise(biomass = mean(biomass), .groups = "drop")




# Function to select sites that have been consistently monitored over the years
select_consistent_sites <- function(data, min_years_threshold = NULL) {
  # Get the total number of years in the dataset
  total_years <- length(unique(data$year))
  message(paste0("Total years in dataset: ", total_years))
  
  # Count how many years each reef has been monitored in each region
  reef_monitoring_consistency <- data %>%
    dplyr::group_by(.data$region, .data$reef) %>%
    dplyr::summarise(
      years_monitored = dplyr::n_distinct(.data$year),
      years_coverage = .data$years_monitored / total_years,
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$region, dplyr::desc(.data$years_monitored))
  
  # If no threshold provided, find the maximum possible coverage
  if (is.null(min_years_threshold)) {
    # Get the maximum coverage for each region
    max_coverage_by_region <- reef_monitoring_consistency %>%
      dplyr::group_by(.data$region) %>%
      dplyr::summarise(max_coverage = max(.data$years_coverage), .groups = "drop")
    
    min_max_coverage <- min(max_coverage_by_region$max_coverage)
    min_years_count <- ceiling(min_max_coverage * total_years)
    
    message(paste0("Maximum consistent coverage across all regions: ", 
                  round(min_max_coverage * 100, 1), "% (", min_years_count, " years)"))
    min_years_threshold <- min_years_count
  }
  
  # Select reefs that have been monitored for at least the threshold number of years
  consistent_reefs <- reef_monitoring_consistency %>%
    dplyr::filter(.data$years_monitored >= min_years_threshold)
  
  # Count how many consistent sites we have per region
  consistent_sites_count <- consistent_reefs %>%
    dplyr::group_by(.data$region) %>%
    dplyr::summarise(site_count = dplyr::n(), .groups = "drop")
  
  message("Number of consistent sites by region (monitored for at least ", 
         min_years_threshold, " years):")
  print(consistent_sites_count)
  
  # Filter data to only include consistently monitored reefs
  filtered_data <- data %>%
    dplyr::inner_join(consistent_reefs %>% dplyr::select(.data$region, .data$reef), 
               by = c("region", "reef"))
  
  return(filtered_data)
}

# Apply the function to get the consistently monitored sites
balanced_biomass <- select_consistent_sites(biomass_raw)

# Verify number of sites per year
sites_check <- balanced_biomass %>%
  dplyr::group_by(.data$year) %>%
  dplyr::summarise(n_sites = dplyr::n_distinct(.data$reef))

print(sites_check)

# Calculate mean biomass by year and region
biomass_by_region <- balanced_biomass %>%
  dplyr::group_by(.data$year, .data$region) %>%
  dplyr::summarise(
    mean_biomass = mean(.data$biomass, na.rm = TRUE),
    sd_biomass=  sd(.data$biomass, na.rm = TRUE),
    sqrt_n=sqrt(dplyr::n()),
    se_biomass = sd(.data$biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop"
  )

# Calculate mean biomass by year and protection level
biomass_by_protection <- balanced_biomass %>%
  dplyr::group_by(.data$year, .data$protection_level) %>%
  dplyr::summarise(
    mean_biomass = mean(.data$biomass, na.rm = TRUE),
    se_biomass = sd(.data$biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop",
    protection_level=factor(protection_level, levels=c("Prohibited", "Allowed", "Open Area"),
                            labels=c("Sin Pesca", "Multi-uso", "Sin Protección"))
  )


## Calculate mean biomass by year and protection level -------

biomass_f <- fish_productivity_data %>%
  group_by(year, region, protection_level, reef, depth, transect, functional_name) %>%
  summarise(biomass = sum(biomass), .groups = "drop") %>%
  group_by(year, region, protection_level, reef,  functional_name) %>%
  summarise(biomass = mean(biomass), .groups = "drop")


# Apply the function to get the consistently monitored sites
balanced_biomass_f <- select_consistent_sites(biomass_f)




biomass_by_fucntional <- balanced_biomass_f %>%
  dplyr::group_by(.data$year, .data$region, .data$protection_level, .data$functional_name) %>%
  dplyr::summarise(
    mean_biomass = mean(.data$biomass, na.rm = TRUE),
    se_biomass = sd(.data$biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop",
    protection_level=factor(protection_level, levels=c("Prohibited", "Allowed", "Open Area"),
                            labels=c("Sin Pesca", "Multi-uso", "Sin Protección"))
  )
### Bar plots -------
abs.biom.2025 <- biomass_by_fucntional %>% 
  filter(year==2025) %>% 
  filter(region!="La Ventana", ) %>% 

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
  filter(functional_name!="Pelágicos", !is.na(functional_name)) %>% 
  group_by(region, functional_name) %>% 
  summarise(biomass=mean(mean_biomass, na.rm=T)) %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=biomass, fill=functional_name, col=functional_name))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
  labs(y="", x="Región", fill="Grupo Funcional", title="2025")+
  theme_bw()+ 

  # coord_flip()+
  theme(axis.text.x = element_text(angle=0),
        axis.title = element_text(face="bold", size=10),
        legend.position = "",
        plot.title = element_text(face="bold", size=14, hjust=0.5),
        legend.title = element_text(hjust=0.5, face="bold"))+
  ylim(0,2.25)



hist.biom.abs <- biomass_by_fucntional %>% 
  # filter(year==2025) %>% 
  filter(!region %in% c("La Ventana","Cabo Pulmo" ) )%>% 
  
  
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
  filter(functional_name!="Pelágicos", !is.na(functional_name)) %>% 
  group_by(region, functional_name) %>% 
  summarise(biomass=mean(mean_biomass, na.rm=T)) %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  ggplot(aes(x=reorder(region, -lat), y=biomass, fill=functional_name, col=functional_name))+
  geom_bar(stat="identity")+
  scale_fill_viridis_d()+ scale_color_viridis_d(guide = "none")+
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



ggsave("figures/biomass_bars_2025_region.png", biomass.bars, width = 12, height = 6, dpi = 300)






#### Bar plot Top 10 Biomass -------------

# Filter data for Asteroidea and Echinoidea
fish.spp <- fish_productivity_data %>%

  # Calculate richness at the reef level
  group_by(year, region, protection_level, reef, depth, transect, functional_name, species) %>%
  
  summarise(biomass = sum(biomass, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, region, protection_level, reef, functional_name, species) %>%
  summarise(biomass = mean(biomass, na.rm = TRUE), .groups = "drop")
# Apply the function to get the consistently monitored sites
balanced_fish <- select_consistent_sites(fish.spp) %>% 
  mutate(period=case_when( year< 2014 ~"Histórico (1998-2013)",
                           year>= 2014 & year<=2024 ~"Calentamiento (2014-2024)",
                           year ==2025 ~"2025"))


functional.fish <- balanced_fish %>% 
  select(species, functional_name) %>% 
  distinct()


 p3 <-  balanced_fish %>% 
    left_join(lats %>% janitor::clean_names()) %>% 
    mutate(degree=ifelse(round(lat,0)>=26, "Arriba de 26°", "Debajo de 26°")) %>%  
    filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
    
  filter(!is.na(functional_name)) %>% 
  # filter(period=="Histórico") %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree!="Debajo de 26°") %>% 
    filter(functional_name!="Pelagic") %>% 

  
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm=T)) %>%
  # filter(MPA=="MPA") %>%
  mutate(period=factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025"))) %>%
   ungroup() %>% 
  group_by(period) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  
   
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
  filter(functional_name!="Pelágicos") %>% 
  ungroup() %>% 
  mutate(species= reorder_within(species, biomass, period)) %>% 
  ggplot(aes(x=species, y = biomass, fill=functional_name)) +
  geom_col()+
    scale_x_continuous(breaks = seq(0,3.5, 0.5))+
  coord_flip()+
  facet_wrap(~period,  scales="free_y", nrow=1)+

  scale_x_reordered() +

  # scale_color_material_d() +
  # scale_fill_manual(values = c("firebrick", "darkgreen"))+
  labs(y = "Biomasa (T/ha)", x="", fill="Grupo Funcional", title="Latitud mayor a 26°") +
  theme_classic() +
  theme(legend.position = "") +
  guides(colour = "none")+
  theme(axis.text.y = element_text(face= "italic"),
        axis.title = element_text(face="bold"),
        plot.title = element_text(hjust=0.5, face="bold", size=15),
        strip.background = element_blank(),  # Removes the background
        strip.text = element_text(face = "bold", hjust = 0.5, size=12),
        legend.title.position = "bottom",
        legend.title = element_text(face="bold", hjust=0.5)
        
  )+
  scale_fill_viridis_d()

  ggsave("figures/top10_fish_biomass_periods.jpg", width = 12, height = 6, dpi = 300)
              

p4 <- 
  balanced_fish %>% 
  left_join(lats %>% janitor::clean_names()) %>% 
  mutate(degree=ifelse(round(lat,0)>=26, "Arriba de 26°", "Debajo de 26°")) %>%  
  filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>% 
  
  filter(!is.na(functional_name)) %>% 
  # filter(period=="Histórico") %>% 
  filter(!region %in% c("Cabo Pulmo", "La Ventana"), degree!="Arriba de 26°") %>% 
  filter(functional_name!="Pelagic") %>% 
  
  
  group_by(period, species) %>% 
  summarise(biomass = mean(biomass, na.rm=T)) %>%
  # filter(MPA=="MPA") %>%
  mutate(period=factor(period, levels = c("Histórico (1998-2013)", "Calentamiento (2014-2024)", "2025"))) %>%
  ungroup() %>% 
  group_by(period) %>%
  top_n(10, biomass) %>%
  left_join(functional.fish) %>% 
  
  
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
  filter(functional_name!="Pelágicos") %>% 
  ungroup() %>% 
  mutate(species= reorder_within(species, biomass, period)) %>% 
  ggplot(aes(x=species, y = biomass, fill=functional_name)) +
  geom_col()+
  scale_x_continuous(breaks = seq(0,3.5, 0.5))+
  coord_flip()+
  facet_wrap(~period,  scales="free_y", nrow=1)+
  
  scale_x_reordered() +
  
  # scale_color_material_d() +
  # scale_fill_manual(values = c("firebrick", "darkgreen"))+
  labs(y = "Biomasa (T/ha)", x="", fill="Grupo Funcional", title="Latitud menor a 26°") +
  theme_classic() +
  theme(legend.position = "bottom") +
  guides(colour = "none")+
  theme(axis.text.y = element_text(face= "italic"),
        axis.title = element_text(face="bold"),
        plot.title = element_text(hjust=0.5, face="bold", size=15),
        strip.background = element_blank(),  # Removes the background
        strip.text = element_text(face = "bold", hjust = 0.5, size=12, color="white"),
        legend.title.position = "bottom",
        legend.title = element_text(face="bold", hjust=0.5))+
  scale_fill_viridis_d()

plot_grid(
  p3 + theme(legend.position = "none"),  # Remove legend from p1
  p4 + theme(legend.position = "bottom"),  # Remove legend from p2
  labels = c("A", "B"),   
  
  # rel_heights = c(3, 0.1),
  # Add labels A and B
  ncol = 1                               # Arrange plots in one column
)




ggsave("figures/top10_fish_biomass_hist.jpg", width = 12, height = 12, dpi = 300)


## Plot trends by region ----
biomass_by_region %>% 
  filter(!region%in% c("Cabo Pulmo", "La Ventana")) %>% 

ggplot( aes(x = year, y = mean_biomass, color = region)) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept =4, linetype=2, col="firebrick",linewidth=1)+
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
        axis.title.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(angle = 0, face="bold"),
        axis.text.x=element_text(angle=90)) +
  scale_y_continuous(labels = comma) +
  scale_color_viridis_d()


ggsave("figures/biomass_trends_by_region_2025_2.png", width = 10, height = 6, dpi = 300)
## Plot trends by protection level ----------
ggplot(biomass_by_protection, aes(x = year, y = mean_biomass, color = protection_level)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass, 
                   ymax = mean_biomass + se_biomass),
               width = 0.2) +
  labs(x = "Año", y = "Biomasa Promedio (T/ha)") +
  scale_x_continuous(breaks = seq(1998, 2025, by = 1)) +
  geom_hline(yintercept =4, linetype=2, col="firebrick",linewidth=1)+
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_text(face="bold"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(angle = 0, face="bold"),
        axis.text.x=element_text(angle=90)
         ) +
  scale_y_continuous(labels = comma) +
  scale_color_viridis_d()

ggsave("figures/biomass_trends_by_protection_2025.png",  width = 10, height = 6, dpi = 300)


## Plot by functional group--------



# Reorder Region factor levels based on latitude
# Invertir para que regiones del norte aparezcan arriba en facet_wrap
region_levels <- rev(lats$Region)

biomass_by_fucntional %>% 
  filter(! region%in% c("Cabo Pulmo", "La Ventana") )%>% 
  # left_join(lats) %>% 
  mutate(region = factor(region, levels = c("Santa Rosalia",
                                            "San Basilio",
                                            "Loreto",
                                            "Corredor",
                                            "La Paz")),
         
 functional_name=factor(functional_name, levels=c("GenPred_solitary",
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
  filter(functional_name!="Pelágicos", !is.na(functional_name)) %>% 
  dplyr::group_by(year, region, functional_name) %>%
  dplyr::summarise(
    mean_biomass = mean(.data$mean_biomass, na.rm = TRUE),
    date = as.Date(paste0(year, "-01-01"))) %>% 
  ggplot(aes(x = date, y = mean_biomass, col = functional_name)) +
  labs(x = "Año", y = "Biomasa (ton/ha)", col="Grupo Funcional") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              method.args = list(family = Gamma(link = "log")),
              fill =NA) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  facet_wrap(~region, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(hjust = 0.5, face="bold"),
        legend.title.position = "top",
        strip.background = element_rect(color=NA, fill=NA),
        strip.text = element_text(face="bold", size=10),
        legend.position = "right") +
  ylim(0, 1
       )+
  
  scale_color_viridis_d()

ggsave("figures/biomass_trends_region_functional_2025.jpg", width = 10, height = 12, dpi = 300)

biomass_by_fucntional %>% 
  # filter(!region %in% c("Cabo Pulmo", "La Ventana")) %>%
  # left_join(lats) %>% 
  mutate(region = factor(region, levels = region_levels),
         
         functional_name=factor(functional_name, levels=c("GenPred_solitary",
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
  filter(functional_name!="Pelágicos", !is.na(functional_name)) %>% 
  dplyr::group_by(year,  functional_name) %>%
  dplyr::summarise(
    mean_biomass = mean(.data$mean_biomass, na.rm = TRUE),
    date = as.Date(paste0(year, "-01-01"))) %>% 
  ggplot(aes(x = date, y = mean_biomass, col = functional_name)) +
  labs(x = "Año", y = "Biomasa (ton/ha)", col="Grupo Funcional") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3),
              method.args = list(family = Gamma(link = "log")),
              fill =NA) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  # facet_wrap(~region, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 12, face="bold"),
        legend.title = element_text(hjust = 0.5, face="bold"),
        legend.title.position = "top",
        strip.text = element_text(face="bold", size=10),
        legend.position = "right") +
  ylim(0, 1
  )+
  scale_color_viridis_d()

ggsave("figures/biomass_trends_functional_2025.jpg",  width = 10, height = 6, dpi = 300)

## Rainplots -----

balanced_biomass_f %>%
  filter(region!="La Ventana") %>% 
 group_by( region, reef,  functional_name) %>%
  # filter(!(region=="Corredor" & biomass >2.5)) %>% 
 summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop",
    protection_level=factor(protection_level, levels=c("Prohibited", "Allowed", "Open Area"),
                            labels=c("Sin Pesca", "Multi-uso", "Sin Protección"))
  ) %>% 
  mutate(
    region = factor(region, levels = region_levels),
    functional_name = factor(
      functional_name,
      levels = rev(c(
        "GenPred_solitary",
        "GenPred_schooling",
        "EpiBent_schooling",
        "Crip_schooling",
        "Crip_solitary",
        "Plank",
        "Pelagic"
      )),
      labels = rev(c(
        "Depredadores solitarios",
        "Depredadores en cardúmenes",
        "Omnívoros en cardúmen",
        "Herbívoros en cardúmen",
        "Crípticos solitarios",
     
        "Planctívoros",
        "Pelágicos"
      ))
    ), 
    log_biom=log1p(mean_biomass)
  ) %>%
  filter(!is.na(functional_name), functional_name != "Pelágicos") %>%
  ggplot(aes(x = functional_name, y = log1p(mean_biomass), fill = functional_name, col=functional_name)) +
  geom_violin(trim = TRUE, alpha=0.5) +
  # geom_jitter(alpha=0.25)+
  geom_boxplot( alpha=0.5, width=0.05)+
  # ylim(0, 3) +
  facet_wrap(~region, nrow = 1) +
  labs(x="Grupo Funcional", y="log(Biomasa + 1)")+
  coord_flip()+
  theme_bw() +
  theme(legend.position = "",
        legend.title = element_blank(),
        axis.title.y = element_text(face="bold"),
        strip.text = element_text(face="bold", size=10),
       
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_text(angle = 0, face="bold"),
        # axis.text.x=element_text(angle=90)
  ) +
  # scale_y_continuous(labels = comma) +
  scale_fill_viridis_d(direction = -1)+  scale_color_viridis_d(direction = -1)
  

ggsave("figures/biomass_rainplots_functional_2025_log.jpg",  width = 12, height = 6, dpi = 300)


# 2025
balanced_biomass_f %>%
  filter(year==2025) %>% 
  filter(region!="La Ventana") %>% 
  group_by( region, reef,  functional_name) %>%
  # filter(!(region=="Corredor" & biomass >2.5)) %>% 
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop",
    protection_level=factor(protection_level, levels=c("Prohibited", "Allowed", "Open Area"),
                            labels=c("Sin Pesca", "Multi-uso", "Sin Protección"))
  ) %>% 
  mutate(
    region = factor(region, levels = region_levels),
    functional_name = factor(
      functional_name,
      levels = rev(c(
        "GenPred_solitary",
        "GenPred_schooling",
        "EpiBent_schooling",
        "Crip_schooling",
        "Crip_solitary",
        "Plank",
        "Pelagic"
      )),
      labels = rev(c(
        "Depredadores solitarios",
        "Depredadores en cardúmenes",
        "Omnívoros en cardúmen",
        "Herbívoros en cardúmen",
        "Crípticos solitarios",
        
        "Planctívoros",
        "Pelágicos"
      ))
    )
  ) %>%
  left_join(lats %>% janitor::clean_names()) %>% 
  filter(!is.na(functional_name), functional_name != "Pelágicos") %>%
  ggplot(aes(x = functional_name, y = log1p(mean_biomass), fill = functional_name, col=functional_name)) +
  geom_violin(trim = TRUE, alpha=0.5) +

  # geom_jitter(alpha=0.25)+
  geom_boxplot( alpha=0.5, width=0.05)+
  # ylim(0, 3) +
  facet_wrap(~region, nrow = 1) +
  labs(x="Grupo Funcional", y="Biomasa Promedio (T/ha)")+
  coord_flip()+
  theme_bw() +
  theme(legend.position = "",
        legend.title = element_blank(),
        axis.title.y = element_text(face="bold"),
        strip.text = element_text(face="bold", size=10),
        
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_text(angle = 0, face="bold"),
        # axis.text.x=element_text(angle=90)
  ) +
  # scale_y_continuous(labels = comma) +
  scale_fill_viridis_d(direction = -1)+  scale_color_viridis_d(direction = -1)


balanced_biomass_f %>%
  filter(region!="La Ventana") %>% 
  group_by( protection_level, reef,  functional_name) %>%
  # filter(!(region=="Corredor" & biomass >2.5)) %>% 
  summarise(
    mean_biomass = mean(biomass, na.rm = TRUE),
    se_biomass = sd(biomass, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups = "drop",
    protection_level=factor(protection_level, levels=c("Prohibited", "Allowed", "Open Area"),
                            labels=c("Sin Pesca", "Multi-uso", "Sin Protección"))
  ) %>% 
  mutate(
    # region = factor(region, levels = region_levels),
    functional_name = factor(
      functional_name,
      levels = rev(c(
        "GenPred_solitary",
        "GenPred_schooling",
        "EpiBent_schooling",
        "Crip_schooling",
        "Crip_solitary",
        "Plank",
        "Pelagic"
      )),
      labels = rev(c(
        "Depredadores solitarios",
        "Depredadores en cardúmenes",
        "Omnívoros en cardúmen",
        "Herbívoros en cardúmen",
        "Crípticos solitarios",
      
        "Planctívoros",
        "Pelágicos"
      ))
    )
  ) %>%
  filter(!is.na(functional_name), functional_name != "Pelágicos") %>%
  ggplot(aes(x = functional_name, y = log1p(mean_biomass), fill = functional_name, col=functional_name)) +
  geom_violin(trim = TRUE, alpha=0.5) +
  # geom_jitter(alpha=0.25)+
  geom_boxplot( alpha=0.5, width=0.05)+
  ylim(0, 3) +
  facet_wrap(~protection_level, nrow = 1) +
  labs(x="Grupo Funcional", y="log(Biomasa + 1)")+
  coord_flip()+
  theme_bw() +
  theme(legend.position = "",
        legend.title = element_blank(),
        axis.title.y = element_text(face="bold"),
        strip.text = element_text(face="bold", size=10),
        
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill=NA, color=NA),
        axis.title.x = element_text(angle = 0, face="bold"),
        # axis.text.x=element_text(angle=90)
  ) +
  # scale_y_continuous(labels = comma) +
  scale_fill_viridis_d(direction = -1)+  scale_color_viridis_d(direction = -1)



ggsave("figures/biomass_rainplots_functional_prot_2025.jpg",  width = 12, height = 6, dpi = 300)





## GAM Fish  --------
# Fit GAM models to assess trends
# By region - using region as factor to avoid 'by' variable error
biomass_by_region$region_factor <- as.factor(biomass_by_region$region)
gam_region <- gam(mean_biomass ~ s(year, by = region_factor) + region_factor, 
                  data = biomass_by_region, 
                  family = gaussian())
summary(gam_region)

# By protection level
biomass_by_protection$protection_factor <- as.factor(biomass_by_protection$protection_level)
gam_protection <- gam(mean_biomass ~ s(year, by = protection_factor) + protection_factor, 
                      data = biomass_by_protection, 
                      family = gaussian())
summary(gam_protection)


## NRSI -------

nrsi <- fish_productivity_data %>% select_consistent_sites()


nrsi  <- fish_productivity_data %>% 

  mutate(index_levels = case_when(
    str_detect(trophic_level_f, "4-4.5") ~ "UTL", 
    str_detect(trophic_level_f, "2-2.5") ~ "LTL", 
    TRUE ~ "CTL"
  ), 
  # reef=ifelse(reef=="CORONADO_PUNTA_LOBOS", "CORONADO_LAJAS", reef)
  # protection = factor(protection, levels = c("Not Protected", "Lightly Protected", "Fish Refuge", "Fully Protected"))
  ) %>%
  mutate(degree = round(latitude)) %>% 
  group_by(year, region,reef,  degree, latitude, longitude, depth2, index_levels, transect) %>% 
  summarise(biomass = sum(biomass, na.rm = TRUE)) %>% 
  group_by(year, region, reef,  index_levels, degree, latitude, longitude) %>% 
  summarise(biomass = mean(biomass, na.rm = TRUE)) %>% 
  group_by(year, region, reef,  degree,) %>% 
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
    mean((sampled$nrsi), na.rm = TRUE)
  })
}


resampled_prod_reef <- nrsi %>%
  
  group_by(region) %>%
  nest() %>%
  mutate(
    ResampledMeans = map(data, ~resample_mean(.x)),
    Median = map_dbl(ResampledMeans, mean)
  ) %>%
  dplyr::select(region, ResampledMeans, Median)



toplot_ltms_reef <- resampled_prod_reef %>%
  unnest(ResampledMeans) %>% 
  left_join(lats %>% janitor::clean_names())



hist.nrsi <- toplot_ltms_reef %>% 
  filter(region!="La Ventana") %>% 


  # filter(sitio %in% c("Coronado Punta Blanca", "Lobera", "Submarino", "La Reinita")) %>% 
  ggplot(aes(y = reorder(region, lat), x = ResampledMeans, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "", y = "", title="Histórico") +
  scale_fill_gradientn(colors = c("firebrick", "orange", "gray", "grey", "green"),
                       values = scales::rescale(c(-1, -0.75, 0, 0.5, 0.75, 1)))+
  xlim(-1, 1)+
  # xlim(-0.5,0.5) +
  theme_bw() +
  theme(legend.position = "", 
        plot.title = element_text(face="bold", size=12, hjust=0.5),
        axis.text.x = element_text(angle = 90, vjust = .5))  


# 2025

resampled_prod_reef <- nrsi %>%
  filter(year==2025) %>% 
  group_by(region) %>%
  nest() %>%
  mutate(
    ResampledMeans = map(data, ~resample_mean(.x)),
    Median = map_dbl(ResampledMeans, mean)
  ) %>%
  dplyr::select(region, ResampledMeans, Median)



toplot_ltms_reef <- resampled_prod_reef %>%
  unnest(ResampledMeans) %>% 
  left_join(lats %>% janitor::clean_names())



nrsi.2025 <- toplot_ltms_reef %>% 
  filter(region!="La Ventana") %>% 
  
  
  # filter(sitio %in% c("Coronado Punta Blanca", "Lobera", "Submarino", "La Reinita")) %>% 
  ggplot(aes(y = reorder(region, lat), x = ResampledMeans, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "NRSI", y = "", title="2025") +
  scale_fill_gradientn(colors = c("firebrick", "orange", "gray", "grey", "grey"),
                       values = scales::rescale(c(-1, -0.75, 0, 0.5, 0.75, 1)))+
  # xlim(-0.5,0.5) +
  theme_bw() +
  xlim(-1, 1)+
  theme(legend.position = "", 
        plot.title = element_text(face="bold", size=12, hjust=0.5),
        axis.text.x = element_text(angle = 90, vjust = .5),
        axis.title = element_text(face="bold"))  



(nrsi.plot <- hist.nrsi + nrsi.2025 + plot_layout(ncol = 1))

ggsave("figures/nrsi_by_region_hist_2025.png",  width = 10, height = 6, dpi = 300)

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

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(writexl)
})

utils::globalVariables(c(
  "IDSpecies","Species","Size","Quantity",
  "Size_min","Size_max","Quantity_min","Quantity_max",
  "n","sample_flag","Type","Label","Year","Month","Day",
  "Region","IDReef","Reef","Depth","Transect"
))

#' Calculate statistical thresholds for outlier detection
#'
#' This function computes robust statistical thresholds for each species
#' using either quantiles (default) or modified Z-scores based on median
#' absolute deviation (MAD), which is more robust to outliers than standard
#' deviations.
#'
#' @param BASE_HISTORICA Dataframe containing monitoring data
#' @param method Character string indicating method: "quantile" or "zscore"
#' @param p Numeric value for quantile probability (default 0.95, i.e., 95% interval)
#' @param z_threshold Threshold for modified z-score (default 3.5, recommended by Iglewicz & Hoaglin)
#'
#' @return Dataframe with thresholds for Size and Quantity per species
outlier_thresholds <- function(BASE_HISTORICA, method = "quantile", p = 0.95, z_threshold = 3.5) {
  # Ensure we have enough observations
  species_counts <- BASE_HISTORICA %>%
    group_by(IDSpecies, Species) %>%
    summarise(n = n(), .groups = "drop")
  
  # Filter for species with sufficient observations
  sufficient_data <- BASE_HISTORICA %>%
    inner_join(filter(species_counts, n >= 5), by = c("IDSpecies", "Species"))
  
  if(nrow(sufficient_data) == 0) {
    warning("No species with sufficient data (n >= 5) found")
    return(NULL)
  }
  
  if(method == "quantile") {
    # Using quantiles (more robust than mean ± z*SD for non-normal distributions)
    scores <- sufficient_data %>%
      group_by(IDSpecies, Species) %>%
      summarise(
        # Sample number
        n = n(),
        
        # Size metrics
        Size_median = median(Size, na.rm = TRUE),
        Size_min = quantile(Size, probs = (1-p)/2, na.rm = TRUE),
        Size_max = quantile(Size, probs = 1-(1-p)/2, na.rm = TRUE),
        Size_IQR = IQR(Size, na.rm = TRUE),
        
        # Quantity metrics 
        Quantity_median = median(Quantity, na.rm = TRUE),
        Quantity_min = max(0, quantile(Quantity, probs = (1-p)/2, na.rm = TRUE)), # Can't be negative
        Quantity_max = quantile(Quantity, probs = 1-(1-p)/2, na.rm = TRUE),
        Quantity_IQR = IQR(Quantity, na.rm = TRUE),
        
        .groups = "drop"
      )
  } else if(method == "zscore") {
    # Using modified Z-scores based on MAD (Median Absolute Deviation)
    # This is particularly robust against outliers
    scores <- sufficient_data %>%
      group_by(IDSpecies, Species) %>%
      summarise(
        # Sample number
        n = n(),
        
        # Size metrics
        Size_median = median(Size, na.rm = TRUE),
        Size_mad = mad(Size, na.rm = TRUE),
        Size_min = Size_median - (z_threshold * Size_mad),
        Size_max = Size_median + (z_threshold * Size_mad),
        Size_IQR = IQR(Size, na.rm = TRUE),
        
        # Quantity metrics
        Quantity_median = median(Quantity, na.rm = TRUE),
        Quantity_mad = mad(Quantity, na.rm = TRUE),
        Quantity_min = max(0, Quantity_median - (z_threshold * Quantity_mad)), # Can't be negative
        Quantity_max = Quantity_median + (z_threshold * Quantity_mad),
        Quantity_IQR = IQR(Quantity, na.rm = TRUE),
        
        .groups = "drop"
      )
  } else {
    stop("Method must be either 'quantile' or 'zscore'")
  }
  
  # Add sample size flags
  refs_meta <- scores %>%
    mutate(
      sample_flag = case_when(
        n < 10 ~ "Very limited data",
        n < 30 ~ "Limited data",
        TRUE ~ "Sufficient data"
      )
    )
  
  return(refs_meta)
}

#' Generate reference thresholds from historical monitoring data
#'
#' This function computes statistical thresholds for each species and saves them to a file
#' for later use in outlier detection.
#'
#' @param HISTORICAL_DB Dataframe containing historical monitoring data
#' @param output_path File path to save the reference thresholds
#' @param method Character string for method: "quantile" (default) or "zscore"
#' @param p Numeric value for quantile probability (default 0.95, i.e., 95% interval)
#'
#' @return Invisibly returns the reference thresholds dataframe
generate_reference_thresholds <- function(HISTORICAL_DB, 
                                      output_path = "data/auxiliar/ref_thresholds.RDS", 
                                      method = "quantile", 
                                      p = 0.95) {
  
  # Calculate robust thresholds
  ref_thresholds <- outlier_thresholds(HISTORICAL_DB, method = method, p = p)
  
  if(is.null(ref_thresholds)) {
    stop("Could not generate thresholds - insufficient data")
  }
  
  # Save for future use
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(ref_thresholds, output_path)
  
  message(paste0("Reference thresholds saved to ", output_path))
  
  # Also save species-specific distribution statistics for visualization
  species_stats <- HISTORICAL_DB %>%
    group_by(IDSpecies, Species) %>%
    summarise(
      size_min = min(Size, na.rm = TRUE),
      size_q25 = quantile(Size, 0.25, na.rm = TRUE),
      size_median = median(Size, na.rm = TRUE),
      size_q75 = quantile(Size, 0.75, na.rm = TRUE),
      size_max = max(Size, na.rm = TRUE),
      quantity_min = min(Quantity, na.rm = TRUE),
      quantity_q25 = quantile(Quantity, 0.25, na.rm = TRUE),
      quantity_median = median(Quantity, na.rm = TRUE),
      quantity_q75 = quantile(Quantity, 0.75, na.rm = TRUE),
      quantity_max = max(Quantity, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    )
  
  # saveRDS(species_stats, gsub("\.RDS$", "_distributions.RDS", output_path))
  
  invisible(ref_thresholds)
}

data_check <- function(new_data, thresholds_path = "data/auxiliar/ref_thresholds.RDS", type = "Both", method = "quantile", historical_data = NULL, p = 0.95, review_dir = dirname(thresholds_path)) {
  if (file.exists(thresholds_path)) {
    ref <- readRDS(thresholds_path)
  } else {
    if (is.null(historical_data)) stop("Reference thresholds not found and historical_data not provided")
    ref <- outlier_thresholds(historical_data, method = method, p = p)
  }
  df <- new_data %>% 
    left_join(ref, by = c("IDSpecies", "Species"))
  size_out <- NULL
  quantity_out <- NULL
  if (type %in% c("Size", "Both") && all(c("Size_min", "Size_max") %in% names(df))) {
    size_out <- df %>%
      filter(!is.na(Size), !is.na(Size_min), (Size < Size_min | Size > Size_max)) %>%
      mutate(Type = "Size") %>%
      select(any_of(c("ID", "Type", "Label", "Year", "Month", "Day", "Region", "IDReef", "Reef", "Depth", "Transect", "IDSpecies", "Species", "Size", "Size_min", "Size_max", "Quantity", "sample_flag")))
  }
  if (type %in% c("Quantity", "Both") && all(c("Quantity_min", "Quantity_max") %in% names(df))) {
    quantity_out <- df %>%
      filter(!is.na(Quantity), !is.na(Quantity_min), (Quantity < Quantity_min | Quantity > Quantity_max)) %>%
      mutate(Type = "Quantity") %>%
      select(any_of(c("ID", "Type", "Label", "Year", "Month", "Day", "Region", "IDReef", "Reef", "Depth", "Transect", "IDSpecies", "Species", "Quantity", "Quantity_min", "Quantity_max", "Size", "sample_flag")))
  }
  outliers <- bind_rows(size_out, quantity_out)
  dir.create(dirname(thresholds_path), showWarnings = FALSE, recursive = TRUE)
  dir.create(review_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(size_out)) {
    writexl::write_xlsx(size_out, file.path(review_dir, "size_check.xlsx"))
  }
  if (!is.null(quantity_out)) {
    writexl::write_xlsx(quantity_out, file.path(review_dir, "quantity_check.xlsx"))
  }
  outliers
}

visualize_thresholds <- function(data, species_list = NULL, thresholds_path = "data/auxiliar/ref_thresholds.RDS", method = "quantile", historical_data = NULL) {
  if (file.exists(thresholds_path)) {
    ref <- readRDS(thresholds_path)
  } else {
    if (is.null(historical_data)) stop("Reference thresholds not found and historical_data not provided")
    ref <- outlier_thresholds(historical_data, method = method)
  }
  if (is.null(species_list)) {
    species_list <- data %>% count(IDSpecies, Species, sort = TRUE) %>% head(6) %>% pull(Species)
  }
  plots <- lapply(species_list, function(sp) {
    d <- data %>% filter(Species == sp)
    r <- ref %>% filter(Species == sp)
    p1 <- ggplot(d, aes(x = Size)) +
      geom_histogram(bins = 30, fill = "grey70") +
      geom_vline(data = r, aes(xintercept = Size_min), col = "red") +
      geom_vline(data = r, aes(xintercept = Size_max), col = "red") +
      labs(title = sp, x = "Size")
    p2 <- ggplot(d, aes(x = Quantity)) +
      geom_histogram(bins = 30, fill = "grey70") +
      geom_vline(data = r, aes(xintercept = Quantity_min), col = "red") +
      geom_vline(data = r, aes(xintercept = Quantity_max), col = "red") +
      labs(x = "Quantity")
    list(size = p1, quantity = p2)
  })
  plots
}

#' Example usage:
#' 
#' # 1. Generate reference thresholds from historical data
#' # generate_reference_thresholds(ltem)
#' 
#' # 2. Check new data for outliers
#' # outliers <- data_check(new_data, type = "Both")
#' 
#' # 3. Visualize distributions and thresholds for key species
#' # visualize_thresholds(ltem, species_list = c("Acanthaster planci", "Diadema mexicanum"))

#' Example usage:
#' 
#' # 1. Generate reference thresholds from historical data
#' # generate_reference_thresholds(ltem)
#' 
#' # 2. Check new data for outliers
#' # outliers <- data_check(new_data, type = "Both")
#' 
#' # 3. Visualize distributions and thresholds for key species
#' # visualize_thresholds(ltem, species_list = c("Acanthaster planci", "Diadema mexicanum"))

#' Example usage:
#' 
#' # 1. Generate reference thresholds from historical data
#' # generate_reference_thresholds(ltem)
#' 
#' # 2. Check new data for outliers
#' # outliers <- data_check(new_data, type = "Both")
#' 
#' # 3. Visualize distributions and thresholds for key species
#' # visualize_thresholds(ltem, species_list = c("Acanthaster planci", "Diadema mexicanum"))


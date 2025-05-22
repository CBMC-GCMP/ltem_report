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

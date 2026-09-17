# EN: Stage 03 — Metadata checks and outlier detection
# Purpose: Generate reference thresholds and flag outliers in Size/Quantity for QA review.
# Inputs: Latest Stage 02 output from outputs/temp/reefs/; historical data from cfg$historical_dir
# Outputs: Thresholds + outliers to outputs/review/meta_check/
# Assumptions: Historical database is representative; thresholds approximate expected ranges.
#
# ES: Etapa 03 — Revisiones de metadatos y detección de valores atípicos
# Propósito: Generar umbrales de referencia y marcar valores atípicos de Tamaño/Cantidad para revisión.
# Entradas: Última salida de la Etapa 02 en outputs/temp/reefs/; datos históricos en cfg$historical_dir
# Salidas: Umbrales + valores atípicos en outputs/review/meta_check/
# Supuestos: La base histórica es representativa; los umbrales aproximan rangos esperados.

# EN: Load required libraries.
# ES: Carga bibliotecas requeridas.
suppressPackageStartupMessages({
  library(writexl)
})

# EN: Helper function to find the most recent file in a directory matching a pattern.
# ES: Función auxiliar para encontrar el archivo más reciente en un directorio que coincide con un patrón.
.latest_file <- function(dir, pattern) {
  # EN: Helper — return the most recent file in 'dir' matching 'pattern'.
  # ES: Auxiliar — devuelve el archivo más reciente en 'dir' que coincide con 'pattern'.
  files <- dir(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

## EN: Stage 03 main function — generates reference thresholds and flags outliers for QA.
## ES: Función principal de la Etapa 03 — genera umbrales de referencia y marca atípicos para control de calidad.
pipeline_meta_check <- function(cfg) {
  # Prefer latest Stage 02 artifact
  # EN: Locate latest cleaned dataset produced by Stage 02.
  # ES: Ubica el conjunto de datos limpio más reciente producido por la Etapa 02.
  stage2_dir <- file.path(cfg$outputs_dir, "temp", "reefs")
  cleaned_path <- .latest_file(stage2_dir, "^stage-02_ltem-format-check_.*\\.(parquet|rds|csv)$")
  if (is.na(cleaned_path)) {
    stop("No Stage 02 artifact found in ", stage2_dir, ". Run Stage 2 first.")
  }
  # EN: Reader utility for multiple formats (xlsx/xls, rds, csv, parquet).
  # ES: Utilidad de lectura para múltiples formatos (xlsx/xls, rds, csv, parquet).
  read_any <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx","xls")) return(readxl::read_xlsx(path))
    if (ext == "rds") return(readRDS(path))
    if (ext == "csv") return(utils::read.csv(path, stringsAsFactors = FALSE))
    if (ext == "parquet") {
      if (requireNamespace("arrow", quietly = TRUE)) return(arrow::read_parquet(path))
      stop("arrow not installed; cannot read parquet: ", path)
    }
    stop("Unsupported extension: ", ext)
  }
  new_ltem <- read_any(cleaned_path)

  # EN: Load the most recent historical dataset; otherwise use configured fallback if present.
  # ES: Carga el conjunto histórico más reciente; en su defecto usa el respaldo configurado si existe.
  hist_path <- .latest_file(cfg$historical_dir, "^ltem_historic_updated_.*\\.RDS$")
  if (is.na(hist_path) && file.exists(cfg$historical_fallback)) {
    hist_path <- cfg$historical_fallback
  }
  if (is.na(hist_path)) stop("No historical data file found")
  historical <- readRDS(hist_path)

  # EN: Optionally merge PFA records if an auxiliary file is present.
  # ES: Opcionalmente une registros PFA si existe un archivo auxiliar.
  pfa_path <- file.path(cfg$historical_dir, "ltem_pfa_2021-2024.RDS")
  if (file.exists(pfa_path)) {
    pfa <- readRDS(pfa_path)
    historical <- rbind(historical, pfa)
  }

  # EN: Source advanced meta-check helpers if available.
  # ES: Carga funciones avanzadas de meta-check si están disponibles.
  if (file.exists("R/meta_check.R")) source("R/meta_check.R")

  # EN: Create meta_check review dir and initialize thresholds artifact path.
  # ES: Crea el directorio de revisión meta_check e inicializa la ruta del artefacto de umbrales.
  review_dir <- file.path(cfg$outputs_dir, "review", "meta_check")
  dir.create(review_dir, showWarnings = FALSE, recursive = TRUE)
  thresholds_path <- if (exists("save_stage", mode = "function")) save_stage(historical, review_dir, 3, "reference-thresholds", fmt = "rds") else {
    p <- file.path(review_dir, sprintf("stage-03_reference-thresholds_%s.rds", format(Sys.time(), "%Y%m%d-%H%M")))
    saveRDS(NULL, p)
    p
  }
  # EN: Generate and write thresholds based on method and 'p' defined in the config.
  # ES: Genera y escribe los umbrales según el método y 'p' definidos en la configuración.
  generate_reference_thresholds(historical, output_path = thresholds_path, method = cfg$threshold_method, p = cfg$threshold_p)
  if (exists("write_sidecar", mode = "function")) write_sidecar(thresholds_path, "RDS of per-species reference thresholds for Size and Quantity, derived from historical data. Use for reproducible outlier detection (Etapa 03: meta_check).")

  # EN: Compute outliers using new data and thresholds; helpers may write temp files but we capture the data frame.
  # ES: Calcula valores atípicos usando datos nuevos y umbrales; los auxiliares pueden escribir temporales pero capturamos el data frame.
  outliers <- data_check(new_ltem, thresholds_path = thresholds_path, type = "Both", method = cfg$threshold_method, p = cfg$threshold_p, review_dir = tempdir())

  # EN: Save outliers as CSV with standardized naming; attach sidecar description if supported.
  # ES: Guarda los atípicos como CSV con nombre estandarizado; adjunta descripción lateral si es compatible.
  outliers_path <- if (exists("save_stage", mode = "function")) save_stage(outliers, review_dir, 3, "outliers", fmt = "csv") else {
    p <- file.path(review_dir, sprintf("stage-03_outliers_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
    utils::write.csv(outliers, p, row.names = FALSE)
    p
  }
  if (exists("write_sidecar", mode = "function")) write_sidecar(outliers_path, "CSV of observations flagged as outliers (Size/Quantity) versus reference thresholds. Use to guide manual QA review.")

  # Bilingual completion message
  nflag <- if (is.data.frame(outliers)) nrow(outliers) else 0L
  if (exists("stage_log_en_es", mode = "function")) {
    stage_log_en_es(
      3,
      "Metadata checks and outlier detection completed",
      "Revisiones de metadatos y detección de valores atípicos completadas",
      review_dir,
      nflag,
      integer(0),
      "Run Stage 4 to update the historical dataset with cleaned new data.",
      "Ejecuta la Etapa 4 para actualizar la base histórica con los datos limpios."
    )
  }

  # EN: Return references to historical data and paths to threshold and outlier artifacts.
  # ES: Devuelve referencias a datos históricos y rutas a artefactos de umbrales y atípicos.
  list(historical = historical, historical_path = hist_path, thresholds_path = thresholds_path, outliers_path = outliers_path)
}


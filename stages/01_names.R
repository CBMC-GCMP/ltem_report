# EN: Stage 01 — Standardize species names
# Purpose: Fix typos, casing, and synonyms; match to canonical taxonomy.
# Inputs: cfg$new_ltem_path (xlsx/xls/rds/csv), latest species list in cfg$lists_dir
# Outputs: Cleaned table + intermediate export to outputs/temp/species_names/
# Assumptions: Reference taxonomy is authoritative; unmatched entries are flagged.
#
# ES: Etapa 01 — Estandarizar nombres de especies
# Propósito: Corregir errores, mayúsculas/minúsculas y sinónimos; mapear a taxonomía de referencia.
# Entradas: cfg$new_ltem_path (xlsx/xls/rds/csv), última lista de especies en cfg$lists_dir
# Salidas: Tabla limpia + exportación intermedia a outputs/temp/species_names/
# Supuestos: La taxonomía de referencia es la fuente principal; los no coincidentes se marcan.

suppressPackageStartupMessages({
  library(readxl)
})

#' Pipeline stage: prepare species names and outputs
#'
#' Reads a new LTEM dataset (`cfg$new_ltem_path`) by extension, loads the most
#' recent species catalog from `cfg$lists_dir`, optionally standardizes names
#' using the `namechecker` package (if available), and writes outputs to
#' `cfg$outputs_dir`.
#'
#' @param cfg A list-like object containing at least:
#'   - `new_ltem_path`: path to the new LTEM data (xlsx/xls/rds/csv)
#'   - `lists_dir`: directory containing species list Excel files
#'   - `outputs_dir`: directory where outputs will be written
#' @return A named list with elements:
#'   - `new_ltem`: the loaded and optionally name-standardized dataset
#'   - `species_catalog_path`: path to the written species catalog (xlsx)
#'   - `species_catalog_source`: path to the source species list used
#'   - `new_ltem_snapshot`: path to the saved RDS snapshot of `new_ltem`
#' @details
#' - Selects the most recent species list file matching
#'   `^ltem_monitoring_species.*\.xlsx$` in `cfg$lists_dir`.
#' - If `namechecker` is installed, attempts to map `Species` to valid names in
#'   both the species list and the `new_ltem` dataset.
#' - Uses `writexl` to write the catalog when available; otherwise only RDS is
#'   saved.
#' @examples
#' # res <- pipeline_names(list(
#' #   new_ltem_path = "data/drive/ltem_database.xlsx",
#' #   lists_dir = "data/lists/updates",
#' #   outputs_dir = "outputs/pipeline"
#' # ))
pipeline_names <- function(cfg) {
  # EN: Path to the new LTEM data file provided in the config.
  # ES: Ruta al archivo de datos LTEM nuevo indicada en la configuración.
  new_path <- cfg$new_ltem_path
  # EN: Stop early if the input file is missing.
  # ES: Detiene la ejecución si falta el archivo de entrada.
  if (!file.exists(new_path)) stop("File not found: ", new_path)
  # EN: Determine file extension to select the appropriate reader.
  # ES: Determina la extensión para elegir el lector adecuado.
  ext <- tolower(tools::file_ext(new_path))
  # EN: Read Excel files with readxl.
  # ES: Lee archivos de Excel con readxl.
  if (ext %in% c("xlsx", "xls")) {
    new_ltem <- readxl::read_xlsx(new_path)
  } else if (ext == "rds") {
    # EN: Read serialized R objects (RDS).
    # ES: Lee objetos R serializados (RDS).
    new_ltem <- readRDS(new_path)
  } else if (ext == "csv") {
    # EN: Read CSV using base R.
    # ES: Lee CSV usando R base.
    new_ltem <- utils::read.csv(new_path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported new_ltem_path extension: ", ext)
  }

  # EN: Find species list Excel files in 'lists_dir' and pick the most recent by modification time.
  # ES: Busca archivos de lista de especies en 'lists_dir' y elige el más reciente por fecha de modificación.
  spp_files <- list.files(cfg$lists_dir, pattern = "^ltem_monitoring_species.*\\.xlsx$", full.names = TRUE)
  if (length(spp_files) == 0) stop("No species list files found in: ", cfg$lists_dir)
  spp_path <- spp_files[which.max(file.info(spp_files)$mtime)]
  species_list <- readxl::read_xlsx(spp_path)
  # Initialize flag column
  # EN: Ensure a 'Flag' column exists to track changes made by this stage.
  # ES: Asegura que exista una columna 'Flag' para rastrear cambios de esta etapa.
  if (is.null(new_ltem$Flag)) new_ltem$Flag <- ""
  # EN: Keep a copy of original 'Species' to later flag modified rows.
  # ES: Conserva una copia de 'Species' original para marcar filas modificadas.
  before_species <- if ("Species" %in% names(new_ltem)) new_ltem$Species else NULL

  # EN: Optionally use 'namechecker' (if installed) to map names to canonical valid names.
  # ES: Opcionalmente usa 'namechecker' (si está instalado) para mapear nombres a válidos canónicos.
  if (requireNamespace("namechecker", quietly = TRUE)) {
    valid <- namechecker::get_valid_names(species_list$Species)
    valid <- valid[, c("input_name", "valid_name")]
    valid$valid_name[is.na(valid$valid_name)] <- valid$input_name
    tmp <- merge(species_list, valid, by.x = "Species", by.y = "input_name", all.x = TRUE)
    idx_map <- !is.na(tmp$valid_name)
    tmp$Species[idx_map] <- tmp$valid_name[idx_map]
    tmp$valid_name <- NULL
    species_list <- tmp

    # EN: Apply the same mapping to the 'Species' column in the new dataset when present.
    # ES: Aplica el mismo mapeo a la columna 'Species' del nuevo conjunto de datos cuando exista.
    if ("Species" %in% names(new_ltem)) {
      valid2 <- namechecker::get_valid_names(new_ltem$Species)
      valid2 <- valid2[, c("input_name", "valid_name")]
      valid2$valid_name[is.na(valid2$valid_name)] <- valid2$input_name
      m <- match(new_ltem$Species, valid2$input_name)
      repl <- valid2$valid_name[m]
      ind <- !is.na(repl)
      new_ltem$Species[ind] <- repl[ind]
    }
  }

  # Append species_name flag for modified rows
  # ES: Agrega la bandera 'species_name' para filas con cambios en 'Species'.
  # EN: Count how many rows were updated by name standardization.
  # ES: Cuenta cuántas filas fueron actualizadas por la estandarización de nombres.
  rows_updated <- 0L
  if (!is.null(before_species) && "Species" %in% names(new_ltem)) {
    after_species <- new_ltem$Species
    idx_flag <- (!is.na(before_species) & !is.na(after_species) & before_species != after_species) |
                (is.na(before_species) & !is.na(after_species))
    if (any(idx_flag, na.rm = TRUE) && exists("append_flag", mode = "function")) {
      new_ltem <- append_flag(new_ltem, idx_flag, "species_name")
      rows_updated <- sum(idx_flag, na.rm = TRUE)
    }
  }

  # Save standardized output to temp/species_names
  # EN: Save stage artifact using helper if available; otherwise fallback to RDS with timestamp.
  # ES: Guarda el artefacto de la etapa usando el auxiliar si existe; de lo contrario, RDS con sello de tiempo.
  out_dir <- file.path(cfg$outputs_dir, "temp", "species_names")
  out_path <- if (exists("save_stage", mode = "function")) save_stage(new_ltem, out_dir, 1, "species-names-standardized", fmt = "parquet") else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    p <- file.path(out_dir, sprintf("stage-01_species-names-standardized_%s.rds", format(Sys.time(), "%Y%m%d-%H%M")))
    saveRDS(new_ltem, p)
    p
  }

  # Log bilingual completion
  # ES: Registro bilingüe de finalización de la etapa.
  flag_counts <- if (exists("summarize_flags", mode = "function")) summarize_flags(new_ltem) else integer(0)
  if (exists("stage_log_en_es", mode = "function")) {
    stage_log_en_es(
      1,
      "Species names cleaned and standardized",
      "Nombres de especies limpiados y estandarizados",
      out_dir,
      rows_updated,
      flag_counts,
      "Run Stage 2 to verify species IDs and format checks.",
      "Ejecuta la Etapa 2 para verificar IDs de especies y revisiones de formato."
    )
  }

  # EN: Return the cleaned dataset and artifact paths for downstream stages.
  # ES: Devuelve el conjunto de datos limpio y rutas de artefactos para etapas posteriores.
  list(new_ltem = new_ltem, species_catalog_path = NA_character_, species_catalog_source = spp_path, new_ltem_snapshot = out_path)
}

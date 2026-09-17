# EN: Stage 02 — Format and ID checks (species and reefs)
# Purpose: Validate and correct Species/IDSpecies and Reef/IDReef; optional size normalization.
# Inputs: Latest Stage 01 output from outputs/temp/species_names/; species and reef lists in cfg$lists_dir
# Outputs: Cleaned table to outputs/temp/reefs/; review CSVs for transect and size checks
# Assumptions: Reference lists are authoritative; changes are tracked in Flag column.
#
# ES: Etapa 02 — Revisiones de formato e IDs (especies y arrecifes)
# Propósito: Validar y corregir Species/IDSpecies y Reef/IDReef; normalización opcional de tamaños.
# Entradas: Última salida de la Etapa 01 en outputs/temp/species_names/; listas de especies y arrecifes en cfg$lists_dir
# Salidas: Tabla limpia en outputs/temp/reefs/; CSVs de revisión para transectos y tamaños
# Supuestos: Las listas de referencia son la fuente principal; los cambios se registran en la columna Flag.

suppressPackageStartupMessages({
  library(readxl)
})

# EN: Helper to get the most recent file matching a pattern in a directory.
# ES: Auxiliar para obtener el archivo más reciente que coincide con un patrón en un directorio.
# Helper: latest file by pattern
.latest_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NA_character_)
  files[which.max(file.info(files)$mtime)]
}

# EN: Stage 02 main function — validates/corrects species and reef IDs/names using lists.
# ES: Función principal de la Etapa 02 — valida/corrige IDs/nombres de especies y arrecifes usando listas.
pipeline_format_check <- function(cfg) {
  # EN: Path to Stage 01 outputs and latest standardized names artifact.
  # ES: Ruta a las salidas de la Etapa 01 y su artefacto más reciente de nombres estandarizados.
  species_dir <- file.path(cfg$outputs_dir, "temp", "species_names")
  cleaned_names_path <- .latest_file(species_dir, "^stage-01_species-names-standardized_.*\\.(parquet|rds|csv)$")
  # EN: Reader that loads file by extension (xlsx/xls, rds, csv, parquet).
  # ES: Lector que carga el archivo según su extensión (xlsx/xls, rds, csv, parquet).
  read_any <- function(path) {
    if (is.na(path) || !file.exists(path)) return(NULL)
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx","xls")) return(readxl::read_xlsx(path))
    if (ext == "rds") return(readRDS(path))
    if (ext == "csv") return(utils::read.csv(path, stringsAsFactors = FALSE))
    if (ext == "parquet") {
      if (requireNamespace("arrow", quietly = TRUE)) return(arrow::read_parquet(path))
      warning("arrow not installed; cannot read parquet: ", path)
      return(NULL)
    }
    stop("Unsupported extension: ", ext)
  }
  new_ltem <- read_any(cleaned_names_path)
  if (is.null(new_ltem)) {
    # EN: Fallback — read the raw new_ltem from cfg$new_ltem_path if Stage 1 artifact is absent.
    # ES: Alternativa — lee el new_ltem crudo desde cfg$new_ltem_path si falta el artefacto de la Etapa 1.
    ext <- tolower(tools::file_ext(cfg$new_ltem_path))
    if (ext %in% c("xlsx", "xls")) {
      new_ltem <- readxl::read_xlsx(cfg$new_ltem_path)
    } else if (ext == "rds") {
      new_ltem <- readRDS(cfg$new_ltem_path)
    } else if (ext %in% c("csv")) {
      new_ltem <- utils::read.csv(cfg$new_ltem_path, stringsAsFactors = FALSE)
    } else {
      stop("Unsupported new_ltem_path extension: ", ext)
    }
  }

  # EN: Load the latest species and reef reference lists from cfg$lists_dir.
  # ES: Carga las listas de referencia más recientes de especies y arrecifes desde cfg$lists_dir.
  spp_path   <- .latest_file(cfg$lists_dir, "^ltem_monitoring_species.*\\.xlsx$")
  reefs_path <- .latest_file(cfg$lists_dir, "^ltem_monitoring_reefs.*\\.xlsx$")
  stopifnot(!is.na(spp_path), !is.na(reefs_path))
  spp   <- readxl::read_xlsx(spp_path)
  reefs <- readxl::read_xlsx(reefs_path)

  # EN: Optional size list path (from config or default within lists_dir). Used to normalize size categories.
  # ES: Ruta opcional de lista de tamaños (desde config o por defecto en lists_dir). Se usa para normalizar tamaños.
  size_list_path <- if (!is.null(cfg$size_list_path)) cfg$size_list_path else file.path(cfg$lists_dir, "ltem_size_list.csv")
  IDSize <- NULL
  if (file.exists(size_list_path)) {
    IDSize <- utils::read.csv(size_list_path, stringsAsFactors = FALSE)
  }

  # EN: Source helper scripts if present (ID corrections and flag utilities).
  # ES: Carga scripts auxiliares si están presentes (correcciones de ID y utilidades de banderas).
  if (file.exists("R/format_check.R")) source("R/format_check.R")
  if (file.exists("R/names_correction.R")) source("R/names_correction.R")
  if (file.exists("R/flags.R")) source("R/flags.R")

  cleaned <- new_ltem
  # EN: Ensure a 'Flag' column exists to record changes made during this stage.
  # ES: Asegura que exista la columna 'Flag' para registrar cambios hechos en esta etapa.
  if (is.null(cleaned$Flag)) cleaned$Flag <- ""
  flag_before <- cleaned$Flag

  # EN: Apply species ID and name corrections using helper functions when available; add flags for changed rows.
  # ES: Aplica correcciones de ID y nombre de especies con auxiliares cuando existan; agrega banderas para filas con cambios.
  if (exists("speciesid", mode = "function")) {
    before <- if ("IDSpecies" %in% names(cleaned)) cleaned$IDSpecies else NULL
    cleaned <- tryCatch(speciesid(cleaned, spp), error = function(e) cleaned)
    if (!is.null(before) && "IDSpecies" %in% names(cleaned) && exists("append_flag", mode = "function")) {
      idx <- (is.na(before) & !is.na(cleaned$IDSpecies)) | (!is.na(before) & !is.na(cleaned$IDSpecies) & before != cleaned$IDSpecies)
      if (any(idx, na.rm = TRUE)) cleaned <- append_flag(cleaned, idx, "species_id")
    }
  }
  if (exists("speciesnames", mode = "function")) {
    before <- if ("Species" %in% names(cleaned)) cleaned$Species else NULL
    cleaned <- tryCatch(speciesnames(cleaned, spp), error = function(e) cleaned)
    if (!is.null(before) && "Species" %in% names(cleaned) && exists("append_flag", mode = "function")) {
      idx <- (is.na(before) & !is.na(cleaned$Species)) | (!is.na(before) & !is.na(cleaned$Species) & before != cleaned$Species)
      if (any(idx, na.rm = TRUE)) cleaned <- append_flag(cleaned, idx, "species_name")
    }
  }

  # EN: Apply reef ID and name corrections similarly; flag rows where values changed.
  # ES: Aplica correcciones de ID y nombre de arrecife de manera similar; marca filas con cambios.
  if (exists("reefsid", mode = "function")) {
    before <- if ("IDReef" %in% names(cleaned)) cleaned$IDReef else NULL
    cleaned <- tryCatch(reefsid(cleaned, reefs), error = function(e) cleaned)
    if (!is.null(before) && "IDReef" %in% names(cleaned) && exists("append_flag", mode = "function")) {
      idx <- (is.na(before) & !is.na(cleaned$IDReef)) | (!is.na(before) & !is.na(cleaned$IDReef) & before != cleaned$IDReef)
      if (any(idx, na.rm = TRUE)) cleaned <- append_flag(cleaned, idx, "reef_id")
    }
  }
  if (exists("reefsname", mode = "function")) {
    before <- if ("Reef" %in% names(cleaned)) cleaned$Reef else NULL
    cleaned <- tryCatch(reefsname(cleaned, reefs), error = function(e) cleaned)
    if (!is.null(before) && "Reef" %in% names(cleaned) && exists("append_flag", mode = "function")) {
      idx <- (is.na(before) & !is.na(cleaned$Reef)) | (!is.na(before) & !is.na(cleaned$Reef) & before != cleaned$Reef)
      if (any(idx, na.rm = TRUE)) cleaned <- append_flag(cleaned, idx, "reef_name")
    }
  }

  # EN: Normalize 'Size' using optional size list mapping; flag rows where Size changed.
  # ES: Normaliza 'Size' usando el mapeo opcional de tamaños; marca filas donde cambió Size.
  size_changed_idx <- NULL
  if (!is.null(IDSize) && exists("sizeid", mode = "function") && all(c("Label", "IDSize") %in% names(cleaned))) {
    before <- if ("Size" %in% names(cleaned)) cleaned$Size else NULL
    cleaned <- tryCatch(sizeid(cleaned, IDSize), error = function(e) cleaned)
    if (!is.null(before) && "Size" %in% names(cleaned) && exists("append_flag", mode = "function")) {
      size_changed_idx <- (is.na(before) & !is.na(cleaned$Size)) | (!is.na(before) & !is.na(cleaned$Size) & before != cleaned$Size)
      if (any(size_changed_idx, na.rm = TRUE)) cleaned <- append_flag(cleaned, size_changed_idx, "formatting")
    }
  }

  # EN: Optional parity check between INV and PEC transects to find coverage discrepancies.
  # ES: Revisión opcional de paridad entre transectos INV y PEC para detectar discrepancias de cobertura.
  transect_flags <- NULL
  if (exists("unique_trnsct", mode = "function") && exists("compare_trnsct", mode = "function") && exists("flag_trnsct", mode = "function")) {
    inv_rows <- if ("Label" %in% names(cleaned)) cleaned[cleaned$Label == "INV", , drop = FALSE] else cleaned
    pec_rows <- if ("Label" %in% names(cleaned)) cleaned[cleaned$Label == "PEC", , drop = FALSE] else cleaned
    inv  <- tryCatch(unique_trnsct(inv_rows), error = function(e) NULL)
    pec  <- tryCatch(unique_trnsct(pec_rows), error = function(e) NULL)
    if (!is.null(inv) && !is.null(pec)) {
      comp <- tryCatch(compare_trnsct(inv, pec), error = function(e) NULL)
      if (!is.null(comp)) transect_flags <- tryCatch(flag_trnsct(comp), error = function(e) NULL)
    }
  }

  # EN: Save the cleaned dataset artifact for Stage 02; use helper if available else write RDS.
  # ES: Guarda el artefacto del conjunto de datos limpio para la Etapa 02; usa auxiliar si existe o escribe RDS.
  out_dir <- file.path(cfg$outputs_dir, "temp", "reefs")
  out_cleaned <- if (exists("save_stage", mode = "function")) save_stage(cleaned, out_dir, 2, "ltem-format-check", fmt = "parquet") else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    p <- file.path(out_dir, sprintf("stage-02_ltem-format-check_%s.rds", format(Sys.time(), "%Y%m%d-%H%M")))
    saveRDS(cleaned, p)
    p
  }

  # EN: Write review CSVs for transect parity flags and size changes, if they exist.
  # ES: Escribe CSVs de revisión para banderas de transectos y cambios de tamaño, si existen.
  if (!is.null(transect_flags)) {
    review_dir_tf <- file.path(cfg$outputs_dir, "review", "transect_check")
    dir.create(review_dir_tf, showWarnings = FALSE, recursive = TRUE)
    tf_path <- if (exists("save_stage", mode = "function")) save_stage(transect_flags, review_dir_tf, 2, "transect-flags", fmt = "csv") else {
      p <- file.path(review_dir_tf, sprintf("stage-02_transect-flags_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
      utils::write.csv(transect_flags, p, row.names = FALSE)
      p
    }
    if (exists("write_sidecar", mode = "function")) write_sidecar(tf_path, "CSV listing transect coverage discrepancies between INV y PEC. Use to verify missing labels per transect.")
  }
  if (!is.null(size_changed_idx) && any(size_changed_idx, na.rm = TRUE)) {
    review_dir_sz <- file.path(cfg$outputs_dir, "review", "size_check")
    dir.create(review_dir_sz, showWarnings = FALSE, recursive = TRUE)
    size_changes <- tryCatch(cleaned[size_changed_idx, , drop = FALSE], error = function(e) NULL)
    if (!is.null(size_changes)) {
      sz_path <- if (exists("save_stage", mode = "function")) save_stage(size_changes, review_dir_sz, 2, "size-changes", fmt = "csv") else {
        p <- file.path(review_dir_sz, sprintf("stage-02_size-changes_%s.csv", format(Sys.time(), "%Y%m%d-%H%M")))
        utils::write.csv(size_changes, p, row.names = FALSE)
        p
      }
      if (exists("write_sidecar", mode = "function")) write_sidecar(sz_path, "CSV of rows where Size was normalized via IDSize mapping. Review to confirm expected canonical sizes.")
    }
  }

  # EN: Log bilingual completion summary including counts of flags; guide to next stage.
  # ES: Registra un resumen bilingüe con conteos de banderas; guía a la siguiente etapa.
  rows_updated <- sum((flag_before %||% "") != (cleaned$Flag %||% ""))
  flag_counts <- if (exists("summarize_flags", mode = "function")) summarize_flags(cleaned) else integer(0)
  if (exists("stage_log_en_es", mode = "function")) {
    stage_log_en_es(
      2,
      "Format and ID checks completed",
      "Revisiones de formato e ID completadas",
      out_dir,
      rows_updated,
      flag_counts,
      "Run Stage 3 to perform metadata checks and outlier detection.",
      "Ejecuta la Etapa 3 para realizar revisiones de metadatos y detección de valores atípicos."
    )
  }

  # EN: Return cleaned dataset and paths to artifacts used downstream.
  # ES: Devuelve el conjunto de datos limpio y rutas a artefactos usados posteriormente.
  list(cleaned_ltem = cleaned, cleaned_ltem_path = out_cleaned, spp_path = spp_path, reefs_path = reefs_path, transect_flags = transect_flags)
}

# EN: Stage 04 — Update historical dataset
# Purpose: Append cleaned new data to historical database with consistent schema and ordering.
# Inputs: Cleaned new data (from Stage 02) and historical dataset from cfg$historical_dir
# Outputs: Updated historical RDS to cfg$historical_dir with standardized filename
# Assumptions: Column union is used; ordering applied for reproducibility.
#
# ES: Etapa 04 — Actualizar base histórica
# Propósito: Anexar los datos nuevos limpios a la base histórica con esquema y orden consistentes.
# Entradas: Datos nuevos limpios (de la Etapa 02) y base histórica desde cfg$historical_dir
# Salidas: RDS histórico actualizado en cfg$historical_dir con nombre estandarizado
# Supuestos: Se usa la unión de columnas; se aplica orden para reproducibilidad.
suppressPackageStartupMessages({
  # library(tidyverse)
})

# EN: Bind cleaned new rows into historical and write an updated artifact.
# ES: Une filas nuevas limpias a la base histórica y escribe un artefacto actualizado.
# EN: No database writes occur here; only file outputs in cfg$outputs_dir/cfg$historical_dir.
# ES: No hay escrituras a base de datos; solo archivos en cfg$outputs_dir/cfg$historical_dir.
pipeline_update <- function(cfg, cleaned_new, historical) {
  # EN: Basic validation — both inputs must be data frames.
  # ES: Validación básica — ambas entradas deben ser data frames.
  stopifnot(is.data.frame(cleaned_new), is.data.frame(historical))

  # Align columns: union of names, then bind
  # EN: Create the union of column names so both data frames share the same schema.
  # ES: Crea la unión de nombres de columnas para que ambos data frames compartan el mismo esquema.
  all_cols <- union(names(historical), names(cleaned_new))
  add_missing <- function(df, cols) {
    # EN: Add any missing columns to the data frame (filled with NA) and reorder.
    # ES: Agrega las columnas faltantes (rellenas con NA) y reordena.
    miss <- setdiff(cols, names(df))
    if (length(miss) > 0) {
      for (m in miss) df[[m]] <- NA
    }
    df[, cols]
  }

  hist_aligned <- add_missing(historical, all_cols)
  new_aligned  <- add_missing(cleaned_new, all_cols)

  # EN: Row-bind historical with new data to form the updated dataset.
  # ES: Une por filas la base histórica con los datos nuevos para formar el conjunto actualizado.
  updated <- rbind(hist_aligned, new_aligned)

  # Arrange for reproducibility
  # EN: If key columns exist, order rows to ensure stable, reproducible output.
  # ES: Si existen columnas clave, ordena filas para asegurar una salida estable y reproducible.
  ord <- intersect(c("Label","Year","Month","Day","Region","Reef","Transect","Depth","IDSpecies","Species"), names(updated))
  if (length(ord) > 0) {
    key <- updated[ord]
    updated <- updated[do.call(order, key), , drop = FALSE]
  }

  # Write output with standardized naming to historical_dir
  # EN: Ensure output directory exists; write using helper if available, else saveRDS fallback with timestamp.
  # ES: Asegura que exista el directorio de salida; escribe con auxiliar si existe, si no usa saveRDS con sello de tiempo.
  dir.create(cfg$historical_dir, showWarnings = FALSE, recursive = TRUE)
  hist_out_path <- if (exists("save_stage", mode = "function")) save_stage(updated, cfg$historical_dir, 4, "historic-updated", fmt = "rds") else {
    p <- file.path(cfg$historical_dir, sprintf("stage-04_historic-updated_%s.rds", format(Sys.time(), "%Y%m%d-%H%M")))
    saveRDS(updated, p)
    p
  }

  # Bilingual completion message
  # EN: Log completion with count of rows appended from the new dataset.
  # ES: Registra la finalización con el conteo de filas anexadas del nuevo conjunto de datos.
  if (exists("stage_log_en_es", mode = "function")) {
    stage_log_en_es(
      4,
      "Historical dataset updated",
      "Base histórica actualizada",
      cfg$historical_dir,
      nrow(new_aligned),
      integer(0),
      "Pipeline complete.",
      "Pipeline completado."
    )
  }

  # EN: Return the updated data frame and the path to the saved artifact.
  # ES: Devuelve el data frame actualizado y la ruta al artefacto guardado.
  list(updated = updated, updated_path = hist_out_path)
}


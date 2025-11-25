# LTEM Report Pipeline (R)

A reproducible data processing and analysis pipeline for the LTEM program. It standardizes new survey data, validates formats and metadata, and (optionally) updates the historical dataset and figures.

## Overview

- **Stages**
  - 01 Names: standardize species names and flags.
  - 02 Format check: validate and clean current survey data.
  - 03 Meta check: compare against historical/reference thresholds.
  - 04 Update: write updated historical dataset (skipped when `dry_run: true`).
- **Runner**: `02-run_pipeline.R` orchestrates all stages using settings in `config.yml`.
- **Changelog**: see `CHANGELOG.md` for notable changes.

## Repository structure

- **Scripts**
  - `00-unify_dbs.R` — helper to unify multiple raw workbooks (optional).
  - `01-update_species_catalog.R` — update/export species catalog (optional).
  - `02-run_pipeline.R` — main pipeline runner.
  - `03-historical_db_append.R` — append cleaned rows to historical dataset.
  - `04-full_trends_analysis.R` — analysis and figure generation.
  - `R/` — reusable helpers (e.g., load/save list utilities, meta checks, flags).
  - `stages/` — stage definitions (`01_names.R`, `02_format_check.R`, `03_meta_check.R`, `04_update.R`).
- **Data & outputs**
  - `data/` — input/reference lists and raw data (ignored by Git).
  - `outputs/` — cleaned data, temp artifacts, historical outputs (ignored).
  - `figures/` — generated figures (ignored).
  - `reports/` — report documents (ignored).

## Requirements

- **R**: 4.x or newer
- **R packages**:
  `yaml`, `tidyverse`, `readxl`, `readr`, `lubridate`, `writexl`, `dplyr`, `tidyr`,
  `ggplot2`, `stringr`, `tidyselect`, `mgcv`, `patchwork`, `lme4`, `scales`,
  `purrr`, `ggridges`, `tidytext`, `googlesheets4`, `namechecker`, `dafishr`

Install from R:
```r
install.packages(c(
  "yaml","tidyverse","readxl","readr","lubridate","writexl",
  "dplyr","tidyr","ggplot2","stringr","tidyselect",
  "mgcv","patchwork","lme4","scales","purrr","ggridges",
  "tidytext","googlesheets4","namechecker"
))
# If 'dafishr' is not on CRAN, install from its source as appropriate
```

## Configuration (`config.yml`)

The pipeline reads a YAML configuration file. Default path is `config.yml` in the repo root; override with env var `LTEM_PIPELINE_CONFIG`.

Example:
```yaml
# LTEM v2 Pipeline configuration
new_ltem_path: "data/raw/2025/OCT-NOV/ltem_OCT-NOV_2025.xlsx"
historical_dir: "outputs/historical"
historical_fallback: "data/ltem_historic_updated_2025-05-13.RDS"
lists_dir: "data/lists/updates/"
outputs_dir: "outputs"
figures_dir: "figures"
threshold_method: "quantile"
threshold_p: 0.95
# When true, Stage 4 (update) is skipped — no writes are performed
# Set to false to apply updates
dry_run: true
```

## How to run

- **RStudio**
  - Open the project, set working directory to the repo root.
  - Ensure packages are installed.
  - Source the runner:
    ```r
    source("02-run_pipeline.R")
    ```
- **Command line (PowerShell)**
  - From the repo root:
    ```powershell
    Rscript .\02-run_pipeline.R
    ```
  - Use a custom config:
    ```powershell
    $env:LTEM_PIPELINE_CONFIG = "C:\\path\\to\\your_config.yml"
    Rscript .\02-run_pipeline.R
    ```

## Outputs

- **Figures** written to `figures/`.
- **Historical dataset** written to `outputs/historical/` when `dry_run: false`.
- Temporary artifacts and review files under `outputs/`.

## Git hygiene

- Large folders (`data/`, `outputs/`, `figures/`, `reports/`, `connect_db/`) are ignored via `.gitignore` to keep the repository lightweight.

---

## Español — Resumen

- **Qué es**: Canalización para estandarizar, validar y (opcionalmente) actualizar la base histórica LTEM.
- **Ejecución**: `02-run_pipeline.R` usa `config.yml`. Ejecuta 4 etapas (Nombres, Formato, Metadatos, Actualización).
- **Config**: Edita `config.yml` (ver ejemplo arriba). Usa `dry_run: true` para no escribir cambios.
- **Requisitos**: R 4.x y paquetes listados (instala con `install.packages(...)`).
- **Salida**: Figuras en `figures/`; histórico en `outputs/historical/` (cuando `dry_run: false`).
- **Git**: Datos y salidas grandes se excluyen del control de versiones con `.gitignore`.

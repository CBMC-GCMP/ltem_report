# Changelog

All notable changes to this project will be documented in this file.

## 2025-11-25

- Added
  - 00-unify_dbs.R: unify multiple LTEM data sources into a single dataset.
  - 01-update_species_catalog.R: update and export species catalog from latest data.
  - 02-run_pipeline.R: orchestrate end-to-end LTEM processing pipeline.
  - 03-historical_db_append.R: append new surveys to historical database.
  - 04-full_trends_analysis.R: run full analysis and figures generation.
  - R/load_latest_ltem_list.R and R/save_ltem_list.R: helpers for loading/saving LTEM lists.
  - R/size_check.R and helpers/ utilities.

- Changed
  - R/meta_check.R: refine checks and outputs.
  - R/resolve_scientific_names.R: improve name standardization and resolution.
  - config.yml: update configuration to support the new pipeline structure.

- Chore
  - .gitignore expanded to ignore data/, outputs/, figures/, reports/, connect_db/, and env/OS temp files.

Notes

- Data, outputs, figures, and reports are intentionally ignored going forward to keep the repo lightweight.
- Previously tracked large files under those paths were not re-added in this change.

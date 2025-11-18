# ============================================================
 # LTEM Program Data Pipeline — Runner (run_pipeline.R)
 # ============================================================
 # EN — What this script is
 # ------------------------------------------------------------
 # This script orchestrates the LTEM data pipeline by running the
 # main stages in order:
 #   0) Optional taxonomy update (pre-stage, if available)
 #   1) Names
 #   2) Format check
 #   3) Meta check
 #   4) Update (skipped when dry_run = TRUE)
 #
 # You control the pipeline via a YAML configuration file. By default,
 # the script looks for a file named "config.yml" in this folder. You
 # can also point to a different file using the environment variable
 # LTEM_PIPELINE_CONFIG.
 #
 # EN — Requirements
 # ------------------------------------------------------------
 # - R (version 4.x or newer)
 # - Package: yaml
 #   If missing, install from R: install.packages("yaml")
 # - A configuration file (config.yml) in this folder or a custom path
 #
 # EN — Quick start (RStudio)
 # ------------------------------------------------------------
 # 1) Open the project/folder in RStudio and set the working directory
 #    to this "pipeline" folder (Session > Set Working Directory > To Source File Location).
 # 2) Ensure the yaml package is installed: install.packages("yaml")
 # 3) Create or verify a config file named "config.yml" in this folder.
 #    Minimal example (save as config.yml):
 #      figures_dir: "figures"
 #      historical_dir: "historical"
 #      dry_run: true   # set to false to enable the update stage
 # 4) Click the "Source" button in RStudio while this file is open,
 #    or run: source("run_pipeline.R")
 #
 # EN — Quick start (Windows PowerShell)
 # ------------------------------------------------------------
 # 1) Open Windows PowerShell.
 # 2) Change directory to the pipeline folder (example):
 #      cd "C:\\Users\\<you>\\CBMC Dropbox\\Eduardo\\R_projects\\ltem-program\\pipeline"
 # 3) Run:
 #      Rscript .\run_pipeline.R
 #
 # EN — Use a custom config file path
 # ------------------------------------------------------------
 # If your config file is not named config.yml or is in a different location,
 # set the environment variable LTEM_PIPELINE_CONFIG before running:
 #   PowerShell example:
 #     $env:LTEM_PIPELINE_CONFIG = "C:\\path\\to\\your_config.yml"
 #     Rscript .\run_pipeline.R
 #
 # EN — Dry run mode
 # ------------------------------------------------------------
 # If config.yml contains dry_run: true, the pipeline will skip the update
 # stage (no writes). Set dry_run: false to apply updates.
 #
 # EN — Troubleshooting
 # ------------------------------------------------------------
 # - Error about missing 'yaml':
 #     install.packages("yaml")
 # - Config not found:
 #     Ensure config.yml exists in this folder or set LTEM_PIPELINE_CONFIG.
 # - Wrong working directory:
 #     Run from the "pipeline" folder so the script can find the "stages" files.
 # - Need to see logs:
 #     The script prints progress messages like "[1/4] Running names stage...".
 #
 # EN — Outputs
 # ------------------------------------------------------------
 # The script ensures that directories declared in your config (e.g.,
 # figures_dir, historical_dir) exist. The update stage applies changes;
 # when dry_run is TRUE, no changes are written.
 #
 # ============================================================
 # ES — ¿Qué es este script?
 # ------------------------------------------------------------
 # Este script coordina la canalización (pipeline) de datos LTEM
 # ejecutando las etapas principales en orden:
 #   0) Actualización de taxonomía (pre-etapa, si está disponible)
 #   1) Nombres
 #   2) Revisión de formato
 #   3) Revisión de metadatos
 #   4) Actualización (omitida cuando dry_run = TRUE)
 #
 # Controlas el pipeline mediante un archivo de configuración YAML.
 # Por defecto busca "config.yml" en esta carpeta. También puedes
 # indicar otra ruta usando la variable de entorno LTEM_PIPELINE_CONFIG.
 #
 # ES — Requisitos
 # ------------------------------------------------------------
 # - R (versión 4.x o más reciente)
 # - Paquete: yaml
 #   Si falta, instala desde R: install.packages("yaml")
 # - Un archivo de configuración (config.yml) en esta carpeta o en otra ruta
 #
 # ES — Inicio rápido (RStudio)
 # ------------------------------------------------------------
 # 1) Abre el proyecto/carpeta en RStudio y establece el directorio de trabajo
 #    en esta carpeta "pipeline" (Session > Set Working Directory > To Source File Location).
 # 2) Asegúrate de tener instalado el paquete yaml: install.packages("yaml")
 # 3) Crea o verifica un archivo "config.yml" en esta carpeta.
 #    Ejemplo mínimo (guardar como config.yml):
 #      figures_dir: "figures"
 #      historical_dir: "historical"
 #      dry_run: true   # cambia a false para habilitar la etapa de actualización
 # 4) Pulsa "Source" en RStudio con este archivo abierto,
 #    o ejecuta: source("run_pipeline.R")
 #
 # ES — Inicio rápido (Windows PowerShell)
 # ------------------------------------------------------------
 # 1) Abre Windows PowerShell.
 # 2) Cambia al directorio de la carpeta pipeline (ejemplo):
 #      cd "C:\\Users\\<tu_usuario>\\CBMC Dropbox\\Eduardo\\R_projects\\ltem-program\\pipeline"
 # 3) Ejecuta:
 #      Rscript .\run_pipeline.R
 #
 # ES — Usar una ruta de configuración personalizada
 # ------------------------------------------------------------
 # Si tu archivo de configuración no se llama config.yml o está en otra ubicación,
 # define la variable de entorno LTEM_PIPELINE_CONFIG antes de ejecutar:
 #   Ejemplo en PowerShell:
 #     $env:LTEM_PIPELINE_CONFIG = "C:\\ruta\\a\\tu_config.yml"
 #     Rscript .\run_pipeline.R
 #
 # ES — Modo de simulación (dry run)
 # ------------------------------------------------------------
 # Si config.yml contiene dry_run: true, el pipeline omitirá la etapa de
 # actualización (no escribe cambios). Usa dry_run: false para aplicar cambios.
 #
 # ES — Resolución de problemas
 # ------------------------------------------------------------
 # - Error por 'yaml' faltante:
 #     install.packages("yaml")
 # - No se encuentra el archivo de configuración:
 #     Verifica que exista config.yml en esta carpeta o define LTEM_PIPELINE_CONFIG.
 # - Directorio de trabajo incorrecto:
 #     Ejecuta desde la carpeta "pipeline" para que el script encuentre los archivos "stages".
 # - Ver registros:
 #     El script imprime mensajes como "[1/4] Running names stage...".
 #
 # ES — Resultados
 # ------------------------------------------------------------
 # El script crea (si no existen) los directorios definidos en tu configuración
 # (por ejemplo, figures_dir, historical_dir). La etapa de actualización aplica
 # cambios; cuando dry_run es TRUE, no se escribe nada.
 # ============================================================

 if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Please install the 'yaml' package: install.packages('yaml')")
}

# load config
# EN: Read the config file path from environment variable 'LTEM_PIPELINE_CONFIG'; if not set, use "config.yml" in this folder.
# ES: Lee la ruta del archivo de configuración desde la variable de entorno 'LTEM_PIPELINE_CONFIG'; si no está definida, usa "config.yml" en esta carpeta.
env_cfg_path <- Sys.getenv("LTEM_PIPELINE_CONFIG", unset = "config.yml")
# EN: Load the YAML configuration file into the list object 'cfg' so all stages can access settings.
# ES: Carga el archivo de configuración YAML en el objeto lista 'cfg' para que todas las etapas usen estos parámetros.
cfg <- yaml::read_yaml(env_cfg_path)

# source helpers
# EN: If the utilities file exists, load helper functions (e.g., directory creation) used by the pipeline.
# ES: Si existe el archivo de utilidades, carga funciones auxiliares (p. ej., creación de carpetas) usadas por el pipeline.
if (file.exists("R/utils_pipeline.R")) source("R/utils_pipeline.R")

# ensure output dirs
# EN: If a function 'make_dirs' is available, let it create any folders specified by the configuration.
# ES: Si la función 'make_dirs' está disponible, permite que cree las carpetas especificadas por la configuración.
if (exists("make_dirs", mode = "function")) make_dirs(cfg)
# EN: Define a small helper to create directories quietly and recursively (no warnings if they already exist).
# ES: Define un pequeño auxiliar para crear directorios en silencio y de forma recursiva (sin avisos si ya existen).
ensure_dir <- function(p) dir.create(p, showWarnings = FALSE, recursive = TRUE)
# EN: Ensure the 'figures' output directory exists as defined in the config.
# ES: Asegura que exista el directorio de salida 'figures' como se define en la configuración.
ensure_dir(cfg$figures_dir)
# EN: Ensure the 'historical' output directory exists as defined in the config.
# ES: Asegura que exista el directorio de salida 'historical' como se define en la configuración.
ensure_dir(cfg$historical_dir)

# pre-stage: taxonomy update (Stage 00)
# EN: If the optional pre-stage script exists, load it to make its functions (e.g., taxonomy update) available.
# ES: Si existe el script opcional de pre-etapa, cárgalo para habilitar sus funciones (p. ej., actualización taxonómica).
if (file.exists("00-update_species_catalog.R")) source("00-update_species_catalog.R")
# EN: If the taxonomy update function exists, run it with the configuration; errors are caught so the pipeline continues.
# ES: Si existe la función de actualización taxonómica, ejecútala con la configuración; los errores se capturan para que el pipeline continúe.
if (exists("pipeline_taxonomy_update", mode = "function")) {
  # EN: Log that the pre-stage (taxonomy update) is starting.
  # ES: Registra que la pre-etapa (actualización taxonómica) está iniciando.
  cat("[0] Running taxonomy update (pre-stage)...\n")
  # EN: Attempt to run the update; 'silent = TRUE' prevents noisy error output and does not stop execution.
  # ES: Intenta ejecutar la actualización; 'silent = TRUE' evita salida ruidosa de errores y no detiene la ejecución.
  try(pipeline_taxonomy_update(cfg), silent = TRUE)
}

# source stages
# EN: Load Stage 01 (Names) definitions so 'pipeline_names' is available.
# ES: Carga las definiciones de la Etapa 01 (Nombres) para habilitar 'pipeline_names'.
source("stages/01_names.R")
# EN: Load Stage 02 (Format check) definitions so 'pipeline_format_check' is available.
# ES: Carga las definiciones de la Etapa 02 (Revisión de formato) para habilitar 'pipeline_format_check'.
source("stages/02_format_check.R")
# EN: Load Stage 03 (Meta-check) definitions so 'pipeline_meta_check' is available.
# ES: Carga las definiciones de la Etapa 03 (Revisión de metadatos) para habilitar 'pipeline_meta_check'.
source("stages/03_meta_check.R")
# EN: Load Stage 04 (Update) definitions so 'pipeline_update' is available.
# ES: Carga las definiciones de la Etapa 04 (Actualización) para habilitar 'pipeline_update'.
source("stages/04_update.R")

# EN: Log that Stage 01 (Names) is starting.
# ES: Registra que inicia la Etapa 01 (Nombres).
cat("[1/4] Running names stage...\n")
# EN: Run the Names stage using 'cfg'; returns results about name resolution.
# ES: Ejecuta la Etapa de Nombres usando 'cfg'; devuelve resultados sobre la resolución de nombres.
names_res <- pipeline_names(cfg)
# EN: Log that Stage 02 (Format check) is starting.
# ES: Registra que inicia la Etapa 02 (Revisión de formato).
cat("[2/4] Running format check stage...\n")
# EN: Run the Format Check stage; returns cleaned current data and validation info.
# ES: Ejecuta la Etapa de Revisión de Formato; devuelve datos actuales limpiados e información de validación.
fmt_res   <- pipeline_format_check(cfg)
# EN: Log that Stage 03 (Meta-check) is starting.
# ES: Registra que inicia la Etapa 03 (Revisión de metadatos).
cat("[3/4] Running meta-check stage...\n")
# EN: Run the Meta-check stage; compares against historical/metadata and returns summaries.
# ES: Ejecuta la Etapa de Revisión de Metadatos; compara con histórico/metadatos y devuelve resúmenes.
meta_res  <- pipeline_meta_check(cfg)

# EN: Prepare a variable to store the result of the Update stage (will remain NULL if skipped).
# ES: Prepara una variable para almacenar el resultado de la Etapa de Actualización (seguirá en NULL si se omite).
updated_res <- NULL
# EN: If 'dry_run' is TRUE in the config, skip writing changes (no database/files updated).
# ES: Si 'dry_run' es TRUE en la configuración, omite escribir cambios (no se actualizan base de datos/archivos).
if (isFALSE(is.null(cfg$dry_run)) && isTRUE(cfg$dry_run)) {
  # EN: Log that the Update stage is skipped because dry-run mode is active.
  # ES: Registra que la Etapa de Actualización se omite porque el modo de simulación está activo.
  cat("[4/4] Update stage skipped (dry_run = TRUE). No database writes.\n")
} else {
  # EN: Log that Stage 04 (Update) is starting.
  # ES: Registra que inicia la Etapa 04 (Actualización).
  cat("[4/4] Running update stage...\n")
  # EN: Run the Update stage using cleaned current data and historical data to apply changes.
  # ES: Ejecuta la Etapa de Actualización usando datos actuales limpiados y datos históricos para aplicar cambios.
  updated_res <- pipeline_update(cfg, fmt_res$cleaned_ltem, meta_res$historical)
}

# EN: Return all stage outputs invisibly so advanced users can capture them when sourcing the script.
# ES: Devuelve de forma invisible los resultados de todas las etapas para que usuarios avanzados puedan capturarlos al ejecutar el script.
invisible(list(names = names_res, format = fmt_res, meta = meta_res, update = updated_res))

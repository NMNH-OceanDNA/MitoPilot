CLI_SETTING_DEFAULTS <- list(min_mean_depth = "10", accept_linear = "0")

#' Add the command-line tables and export columns that are missing (idempotent)
#' @noRd
migrate_cli_schema <- function(con) {
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS project_settings (key TEXT PRIMARY KEY, value TEXT)")
  for (k in names(CLI_SETTING_DEFAULTS)) {
    DBI::dbExecute(con, "INSERT OR IGNORE INTO project_settings (key, value) VALUES (?, ?)",
      params = list(k, CLI_SETTING_DEFAULTS[[k]]))
  }
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS review_log (
    ID TEXT NOT NULL, path INTEGER, scaffold INTEGER,
    stage TEXT NOT NULL, action TEXT NOT NULL, rules TEXT, reason TEXT, time_stamp INTEGER)")
  cols <- c(start_aa = "INTEGER", stop_aa = "INTEGER", ident_pct = "REAL",
            aa_alignments = "INTEGER", gene_export = "INTEGER")
  defs <- c(start_aa = 10, stop_aa = 10, ident_pct = 60, aa_alignments = 1, gene_export = 0)
  have <- DBI::dbListFields(con, "export_opts")
  for (k in names(cols)) {
    if (!k %in% have) {
      DBI::dbExecute(con, sprintf("ALTER TABLE export_opts ADD COLUMN %s %s", k, cols[[k]]))
      DBI::dbExecute(con, sprintf("UPDATE export_opts SET %s = %s", k, defs[[k]]))
    }
  }
  invisible(NULL)
}

#' Open a project for the command line
#'
#' @return list(path, db, con, dir_out, user_asmb). The caller disconnects `con`.
#' @noRd
open_project <- function(path = ".", migrate = TRUE) {
  path <- normalizePath(path, mustWork = FALSE)
  db <- file.path(path, ".sqlite")
  cfg <- file.path(path, ".config")
  if (!file.exists(db) || !file.exists(cfg)) {
    cli_fail(paste0("Not a MitoPilot project (no .sqlite and .config): ", path))
  }
  conf <- readLines(cfg)
  dir_out <- conf |>
    stringr::str_extract("publishDir.*") |>
    na.omit() |>
    stringr::str_remove("^[^'|^\"]+['\"]") |>
    stringr::str_extract("^[^'|^\"]+")
  if (length(dir_out) == 0) cli_fail(paste0("No publishDir in .config: ", cfg))
  asmb_line <- grep("asmbDir", conf, value = TRUE)[1]
  asmb <- if (is.na(asmb_line)) NA_character_ else stringr::str_trim(stringr::str_split(
    stringr::str_split(asmb_line, "=")[[1]][2], "'")[[1]][2])
  if (migrate) {
    gaps <- local({
      con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db)
      on.exit(DBI::dbDisconnect(con))
      schema_gaps(con)
    })
    if (length(gaps) > 0) backwards_compatibility(path = path, update_config = FALSE)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db)
  options(MitoPilot.db = db)
  list(path = path, db = db, con = con, dir_out = file.path(path, dir_out[1]),
       user_asmb = length(asmb) == 1 && !is.na(asmb) && asmb != "NA")
}

#' @noRd
project_settings <- function(con) {
  s <- DBI::dbGetQuery(con, "SELECT key, value FROM project_settings")
  v <- stats::setNames(as.list(s$value), s$key)
  list(min_mean_depth = as.numeric(v$min_mean_depth),
       accept_linear = identical(as.character(v$accept_linear), "1"))
}

#' @noRd
set_project_setting <- function(con, key, value) {
  if (!key %in% names(CLI_SETTING_DEFAULTS)) {
    cli_fail(paste0("Unknown project setting: ", key))
  }
  if (is.logical(value)) value <- as.integer(value)
  DBI::dbExecute(con, "INSERT OR REPLACE INTO project_settings (key, value) VALUES (?, ?)",
    params = list(key, as.character(value)))
  invisible(NULL)
}

#' Base name of the live run in a project, or NA
#' @noRd
run_active <- function(path) {
  logs <- list.files(run_dir(path, "nextflow"), pattern = "\\.nextflow\\.log$", full.names = TRUE)
  if (length(logs) == 0) return(NA_character_)
  log <- logs[order(run_base_time(run_log_base(logs)), logs, decreasing = TRUE)][1]
  base <- run_log_base(log)
  if (file.exists(file.path(run_dir(path, "reports"), paste0(base, "_report.dcf")))) return(NA_character_)
  if (grepl("No log update for", run_progress_text(log), fixed = TRUE)) return(NA_character_)
  base
}

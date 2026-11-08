#' Create the project directory, normalize input paths, copy the mapping file,
#' and initialize an RStudio project
#'
#' @return list(path, data_path, assembly_path, mapping_out, open_rproj)
#' @noRd
.setup_project_dir <- function(path, mapping_fn, data_path, assembly_path, no_raw_data, Rproj) {
  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  path <- normalizePath(path)
  if (!no_raw_data && length(data_path) == 1) {
    data_path <- normalizePath(data_path, mustWork = FALSE)
  }
  if (length(assembly_path) == 1 && !identical(assembly_path, "NA")) {
    assembly_path <- normalizePath(assembly_path, mustWork = FALSE)
  }
  mapping_out <- file.path(path, "mapping.csv")
  if (!identical(normalizePath(mapping_fn), normalizePath(mapping_out, mustWork = FALSE))) {
    file.copy(mapping_fn, mapping_out, overwrite = TRUE)
  }

  # (optional & only if running from RStudio)
  open_rproj <- FALSE
  if (Rproj && !isFALSE(Sys.getenv("RSTUDIO", FALSE))) {
    if (isFALSE(requireNamespace("rstudioapi", quietly = TRUE))) {
      message("package 'rstudioapi' not available. Skipping RStudio project initialization.")
    } else {
      rstudioapi::initializeProject(path)
      open_rproj <- TRUE
    }
  }
  list(path = path, data_path = data_path, assembly_path = assembly_path,
       mapping_out = mapping_out, open_rproj = open_rproj)
}

#' Resolve a saved profile / built-in template (or use an explicit path),
#' then write it to `<path>/.config` with the per-project placeholders filled
#' @noRd
.write_project_config <- function(path, config, executor, profile_dir, container,
                                  data_path, assembly_path, ncbi_api_key) {
  config <- config %||% resolve_config(executor, profile_dir = profile_dir)
  readLines(config) |>
    fill_config(list(
      CONTAINER_ENGINE = container_engine_block(
        if (executor %in% c("local", "awsbatch")) "docker" else "singularity"),
      CONTAINER_ID = container,
      RAW_DIR = data_path,
      ASMB_DIR = assembly_path,
      MIN_DEPTH = format(2000000, scientific = FALSE),
      NCBI_API_KEY = ncbi_api_key %||% ""
    )) |>
    writeLines(file.path(path, ".config"))
}

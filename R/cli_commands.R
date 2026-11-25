#' @noRd
cli_parse_cmd <- function(args, cmd, desc, flags = character()) {
  p <- cli_parse_args(args, flags = c(flags, "help"))
  if (isTRUE(p$opts$help)) {
    message(paste0(
      "Usage: mitopilot ", cmd, " [options]\n\nOptions:\n",
      paste0("  --", names(desc), "  ", desc, collapse = "\n")
    ))
    return(NULL)
  }
  cli_check_opts(p$opts, c(names(desc), "help"))
  p
}

#' @noRd
cli_num <- function(x, name) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) cli_fail(paste0("--", name, " must be a number, got '", x, "'"))
  v
}

#' @noRd
cli_cmd_init <- function(args) {
  desc <- c(
    test = "make a test project: reads, assemblies, or annotations",
    mapping = "mapping file (CSV)",
    data = "directory of raw reads",
    assemblies = "directory of user assemblies",
    annotations = "directory of user annotations",
    n = "number of test samples",
    "full-size" = "use full-size test data",
    "no-raw-data" = "project has no raw reads",
    executor = "Nextflow executor (local, slurm, ...)",
    container = "container image",
    force = "overwrite an existing project",
    "no-metadata" = "skip online metadata fetch",
    "min-depth" = "minimum read depth",
    "genetic-code" = "NCBI genetic code table",
    settings = "reserved"
  )
  p <- cli_parse_cmd(args, "init", desc,
    flags = c("full-size", "no-raw-data", "force", "no-metadata"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (!is.null(o[["settings"]])) cli_fail("--settings arrives with 'mitopilot opts'")
  if (length(p$pos) > 1L) cli_fail("init takes at most one directory")

  if (!is.null(o[["test"]])) {
    fn <- switch(o[["test"]],
      reads = "new_test_project",
      assemblies = "new_test_project_userAsmb",
      annotations = "new_test_project_userAnnot",
      cli_fail("--test must be one of: reads, assemblies, annotations")
    )
  } else if (!is.null(o[["annotations"]])) {
    fn <- "new_project_userAnnot"
  } else if (!is.null(o[["assemblies"]])) {
    fn <- "new_project_userAsmb"
  } else {
    fn <- "new_project"
    if (is.null(o[["mapping"]]) || is.null(o[["data"]])) cli_fail("init needs --mapping and --data (or --test)")
  }
  fun <- get(fn, envir = asNamespace("MitoPilot"), mode = "function")
  ok <- setdiff(names(formals(fun)), "...")

  a <- list()
  if (length(p$pos)) a$path <- p$pos[1]
  a$mapping_fn <- o[["mapping"]]
  a$data_path <- o[["data"]]
  a$assembly_path <- o[["assemblies"]]
  a$annotation_path <- o[["annotations"]]
  a$n <- if (!is.null(o[["n"]])) cli_num(o[["n"]], "n")
  a$full_size <- if (isTRUE(o[["full_size"]])) TRUE
  a$no_raw_data <- if (isTRUE(o[["no_raw_data"]])) TRUE
  a$executor <- o[["executor"]]
  a$container <- o[["container"]]
  a$force <- if (isTRUE(o[["force"]])) TRUE
  a$min_depth <- if (!is.null(o[["min_depth"]])) cli_num(o[["min_depth"]], "min-depth")
  a$genetic_code <- if (!is.null(o[["genetic_code"]])) cli_num(o[["genetic_code"]], "genetic-code")
  if (isTRUE(o[["no_metadata"]])) {
    if ("fetch_metadata" %in% ok) {
      a$fetch_metadata <- FALSE
    } else {
      a[c("fetch_geome", "fetch_gbif", "fetch_ncbi")] <- list(FALSE, FALSE, FALSE)
    }
  }
  a$Rproj <- FALSE
  bad <- setdiff(names(a), ok)
  if (length(bad)) {
    cli_fail(paste0("'", fn, "' does not accept: ", paste0("--", gsub("_", "-", bad), collapse = ", ")))
  }
  do.call(fun, a)
  0L
}

#' @noRd
cli_cmd_status <- function(args) {
  desc <- c(
    project = "project directory (default: current directory)",
    json = "print JSON instead of a table",
    stage = "show only this stage (assemble, annotate, export)",
    parked = "show only units waiting for a decision"
  )
  p <- cli_parse_cmd(args, "status", desc, flags = c("json", "parked"))
  if (is.null(p)) return(0L)
  o <- p$opts
  proj <- open_project(o$project %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  st <- project_status(proj$con)
  parked <- !st$locked & nzchar(st$reasons)
  show <- st
  keep <- rep(TRUE, nrow(st))
  if (!is.null(o$stage)) keep <- keep & st$stage == o$stage
  if (isTRUE(o$parked)) keep <- keep & parked
  show <- st[keep, , drop = FALSE]
  tally <- warning_tally(proj$con)
  if (isTRUE(o$json)) {
    cat(jsonlite::toJSON(show, dataframe = "rows", na = "null"), "\n", sep = "")
  } else {
    if (nrow(show)) print(show, row.names = FALSE) else cat("No matching units\n")
    message(sprintf("%d units: %d locked, %d parked", nrow(st), sum(st$locked), sum(parked)))
    if (nrow(tally)) {
      message("Warnings: ", paste0(tally$warning, " (", tally$n, ")", collapse = ", "))
    }
  }
  if (any(parked)) 2L else 0L
}

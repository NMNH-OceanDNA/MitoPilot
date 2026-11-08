#' Initialize new MitoPilot Project with user-provided assemblies and annotations
#'
#' Each sample brings one single-sequence mitogenome FASTA and one annotation
#' file (GFF3, BED, or 5-column feature table). The annotation tools are not
#' run: the pipeline curates, validates, and exports the imported features.
#' Gene names found in the files are guessed at setup and must be reviewed in
#' the app before the pipeline runs.
#'
#' @inheritParams new_project_userAsmb
#' @param mapping_fn Path to a mapping file. Should be a csv that minimally
#'   includes an `ID` column with a unique identifier for each sample, a `Taxon`
#'   column containing taxonomic information for each sample, and columns
#'   `R1` and `R2` specifying the names of the raw paired read inputs (unless
#'   `no_raw_data = TRUE`), an `Assembly` column containing names of assembly
#'   fasta files (one sequence each), and an `Annotation` column containing names
#'   of annotation files (`.gff`, `.gff3`, `.bed`, or `.tbl`). An optional
#'   `Topology` column ("circular" or "linear", blank = linear) declares each
#'   assembly's topology; a feature spanning the origin needs "circular". May
#'   include additional columns with other sample metadata.
#' @param annotation_path Path to the directory where the annotation files are
#'   located. Each file must name one seqid equal to the first word of its
#'   FASTA header. Files are copied to `user_annotations/raw/` in the project.
#' @param ... Additional arguments passed as default processing parameters to
#'   [new_db_userAsmb()] through [new_db_userAnnot()]. Mitogenome search,
#'   circularization, and scaffold joining are always off.
#'
#' @export
#'
new_project_userAnnot <- function(
    path = ".",
    mapping_fn = NULL,
    mapping_id = "ID",
    mapping_geome = "GEOME_BCID",
    fetch_geome = TRUE,
    mapping_gbif = "GBIF_ID",
    fetch_gbif = TRUE,
    mapping_ncbi = "BioSample",
    fetch_ncbi = TRUE,
    link_sources = FALSE,
    data_path = NULL,
    no_raw_data = FALSE,
    assembly_path = "NA",
    annotation_path = "NA",
    genetic_code = NULL,
    executor = c("local", "awsbatch", "slurm", "sge", "pbs", "lsf", "NMNH_Hydra", "NOAA_SEDNA"),
    container = paste0("macguigand/mitopilot:", utils::packageVersion("MitoPilot")),
    config = NULL,
    profile_dir = mitopilot_config_dir(),
    ncbi_api_key = NULL,
    Rproj = TRUE,
    force = FALSE,
    ...) {

  check_nextflow_version("new_project_userAnnot")

  executor <- executor[1]
  dots <- list(...)
  forced <- intersect(names(dots), c("find_mitogenome", "mitofinder_db",
                                     "attempt_circularization", "join_scaffolds"))
  if (length(forced) > 0L) {
    stop("Not available for a user-annotation project (always off): ",
         paste(forced, collapse = ", "), call. = FALSE)
  }
  dots$mapping_geome <- mapping_geome
  dots$fetch_geome <- fetch_geome
  dots$mapping_gbif <- mapping_gbif
  dots$fetch_gbif <- fetch_gbif
  dots$mapping_ncbi <- mapping_ncbi
  dots$fetch_ncbi <- fetch_ncbi
  dots$link_sources <- link_sources
  if (no_raw_data) {
    data_path <- "NA"
    message("no_raw_data = TRUE: skipping read mapping and coverage calculation.")
  }
  preflight_project(
    path = path, mapping_fn = mapping_fn, mapping_id = mapping_id,
    data_path = data_path, no_raw_data = no_raw_data, user_asmb = TRUE,
    assembly_path = assembly_path, user_annot = TRUE, annotation_path = annotation_path,
    executor = executor, config = config, profile_dir = profile_dir,
    container = container, genetic_code = genetic_code,
    ncbi_api_key = ncbi_api_key, force = force, db_fun = new_db_userAsmb, dots = dots
  )

  # Project directory, mapping copy, RStudio project ----
  s <- .setup_project_dir(path, mapping_fn, data_path, assembly_path, no_raw_data, Rproj)
  path <- s$path
  data_path <- s$data_path
  assembly_path <- s$assembly_path
  mapping_out <- s$mapping_out
  if (s$open_rproj) on.exit(rstudioapi::openProject(path, newSession = TRUE))
  if (length(annotation_path) == 1 && !identical(annotation_path, "NA")) {
    annotation_path <- normalizePath(annotation_path, mustWork = FALSE)
  }

  # Initialize sqlite db ----
  db <- file.path(path, ".sqlite")
  if (file.exists(db)) {
    message("Overwriting existing database")
    file.remove(db)
  }
  if (isTRUE(force)) {
    unlink(list.files(file.path(path, "user_annotations"), "_annotations\\.csv$",
                      full.names = TRUE))
  }

  imp <- new_db_userAnnot(
    db_path = db,
    mapping_fn = mapping_out,
    mapping_id = mapping_id,
    assembly_path = assembly_path,
    annotation_path = annotation_path,
    no_raw_data = no_raw_data,
    genetic_code = genetic_code,
    mapping_geome = mapping_geome,
    fetch_geome = fetch_geome,
    mapping_gbif = mapping_gbif,
    fetch_gbif = fetch_gbif,
    mapping_ncbi = mapping_ncbi,
    fetch_ncbi = fetch_ncbi,
    link_sources = link_sources,
    ...
  )

  # Original annotation files ----
  mapping <- utils::read.csv(mapping_out)
  raw_dir <- file.path(path, "user_annotations", "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  src <- .ua_path(annotation_path, mapping$Annotation)
  file.copy(src, file.path(raw_dir, paste0(mapping[[mapping_id]], ".",
                                           tolower(tools::file_ext(src)))),
            overwrite = TRUE)

  .write_project_config(path, config, executor, profile_dir, container,
                        data_path, assembly_path, ncbi_api_key)

  message("Project initialized: ", path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  user_annot_summary(con, imp)
  DBI::dbDisconnect(con)
  message("Open the app to review gene names before running the pipeline.")
  message(open_app_hint(path))
}

#' Initialize new MitoPilot Project with user-provided assemblies and annotations
#'
#' Each sample brings one single-sequence mitogenome FASTA and one annotation
#' file (GFF3, BED, or 5-column feature table). The annotation tools are not
#' run: the pipeline curates, validates, and exports the imported features.
#' Gene names found in the files are guessed at setup and must be reviewed in
#' the app before the pipeline runs.
#'
#' @param path Path to the project directory (default = current working
#'   directory)
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
#' @param mapping_id The name of the column in the mapping file that contains
#'   the unique sample identifiers (default = "ID").
#' @param mapping_geome Name of the mapping-file column holding GEOME BCIDs
#'   (optional). Stored as `GEOME_BCID`. See `vignette("Specimen-Metadata")`.
#'   Passed to `new_db_userAnnot()`.
#' @param fetch_geome Fetch GEOME metadata for samples with a BCID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_geome()] later.
#'   Passed to `new_db_userAnnot()`.
#' @param mapping_gbif Name of the mapping-file column holding GBIF occurrence
#'   IDs (optional). Stored as `GBIF_ID`. See `vignette("Specimen-Metadata")`.
#'   Passed to `new_db_userAnnot()`.
#' @param fetch_gbif Fetch GBIF metadata for samples with a GBIF ID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_gbif()] later.
#' @param mapping_ncbi Name of the mapping-file column holding NCBI
#'   BioSample or SRA accessions (optional). Must not be the sample ID column.
#'   Stored as `BioSample`. See `vignette("Specimen-Metadata")`.
#' @param fetch_ncbi Fetch NCBI metadata for samples with a BioSample
#'   value during setup (default TRUE). Set FALSE when offline and run
#'   [fetch_ncbi()] later.
#' @param link_sources Follow links between GEOME, GBIF, and NCBI records to
#'   fill in IDs a sample does not have yet (default FALSE). Saved as the project
#'   setting. See `vignette("Specimen-Metadata")`.
#'   Passed to `new_db_userAnnot()`.
#' @param data_path Path to the directory where the raw data is located. Can be
#'   a AWS s3 bucket even if not using AWS for pipeline execution. Not required
#'   when `no_raw_data = TRUE`.
#' @param no_raw_data (logical) Run without raw sequence data (default = FALSE).
#'   When TRUE, the read-mapping coverage step is skipped: `data_path` is
#'   ignored, coverage/depth statistics are left empty, and annotation coverage
#'   trimming is disabled. Use this to annotate an assembly you already have.
#' @param assembly_path Path to the directory where the mitogenome assemblies are located.
#' @param annotation_path Path to the directory where the annotation files are
#'   located. Each file must name one seqid equal to the first word of its
#'   FASTA header. Files are copied to `user_annotations/raw/` in the project.
#' @param genetic_code Optional NCBI translation table override. Default `NULL`
#'   auto-selects from each sample's curation ruleset; a number sets a
#'   project-wide override. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#' @param executor The executor to use for running the nextflow pipeline. May be
#'   a built-in template ("local" (default), "awsbatch", "slurm", "sge", "pbs",
#'   "lsf", "NMNH_Hydra", "NOAA_SEDNA") or the name of a saved cluster profile
#'   created with [generate_config()]. See [list_configs()] for available names.
#' @param Rproj (logical) Initialize and open an RStudio project in the project
#'   directory (default = TRUE). This option has no effect if not running
#'   interactively in RStudio.
#' @param force (logical) Force recreating of existing project database and
#'   config files (default = FALSE).
#' @param config (optional) provide a path to an existing custom nextflow config
#'   file. If not provided a config file template will be created based on the
#'   specified executor.
#' @param profile_dir Directory searched for saved cluster profiles when
#'   resolving `executor` (default [mitopilot_config_dir()]).
#' @param container The docker container to use for pipeline execution.
#' @param ncbi_api_key Optional NCBI API key string. Used to raise NCBI request
#'   rate limits for the remote BLAST + GenBank fetch steps. See
#'   <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/>. May be left
#'   empty and edited later in `.config` (`params.ncbi_api_key`).
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

  # Create directory if it doesn't exist ----
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
  if (length(annotation_path) == 1 && !identical(annotation_path, "NA")) {
    annotation_path <- normalizePath(annotation_path, mustWork = FALSE)
  }
  mapping_out <- file.path(path, "mapping.csv")
  if (!identical(normalizePath(mapping_fn), normalizePath(mapping_out, mustWork = FALSE))) {
    file.copy(mapping_fn, mapping_out, overwrite = TRUE)
  }

  # Initialize RStudio Project ----
  if (Rproj && !isFALSE(Sys.getenv("RSTUDIO", FALSE))) {
    if (isFALSE(requireNamespace("rstudioapi", quietly = TRUE))) {
      message("package 'rstudioapi' not available. Skipping RStudio project initialization.")
    } else {
      rstudioapi::initializeProject(path)
      on.exit(rstudioapi::openProject(path, newSession = TRUE))
    }
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

  # Config file ----
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

  message("Project initialized: ", path)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  user_annot_summary(con, imp)
  DBI::dbDisconnect(con)
  message("Open the app to review gene names before running the pipeline.")
  message(open_app_hint(path))
}

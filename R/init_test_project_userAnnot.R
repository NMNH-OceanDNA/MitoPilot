#' Initialize a test project for user-supplied annotations
#'
#' Sets up a test project for the user-annotation workflow, with five samples
#' built from the two single-contig assemblies of [new_test_project_userAsmb()]:
#' a linear Conger oceanicus mitogenome with an NCBI-style GFF3 (`UA_LINEAR`)
#' and as BED (`UA_LINEAR_BED`), and a circular Xyrichtys novacula mitogenome
#' with a MITOS-style GFF (`UA_CIRCULAR`), a deliberately messy GFF3
#' (`UA_CIRCULAR_MESSY`), and a MitoPilot export feature table
#' (`UA_CIRCULAR_TBL`). Between them the files cover partial ends, a feature
#' spanning the origin, joined CDS parts, a pseudogene, and names that cannot be
#' mapped. The assemblies and annotations ship with the package; the raw reads
#' are the same fish data used by [new_test_project()].
#'
#' @param path path to the directory for the test project (default = current
#'   working directory). Will be created if it does not already exist.
#' @param no_raw_data (logical) Set the project up without raw reads (default =
#'   FALSE). See [new_project_userAnnot()].
#' @param full_size (logical) Use the full size raw reads (default = FALSE).
#'   Setting to TRUE will download the reads from ENA, which will require
#'   several GB and will take some time to complete. By default the smaller
#'   pre-filtered read files packaged with MitoPilot are used.
#' @param executor The executor to use for running the nextflow pipeline. A
#'   built-in template ("local" (default), "awsbatch", "slurm", "sge", "pbs",
#'   "lsf", "NMNH_Hydra", "NOAA_SEDNA") or a saved profile from
#'   [generate_config()]. See [list_configs()].
#' @param container The container to use for running the pipeline.
#' @param Rproj (logical) Initialize and open an RStudio project in the project
#'   directory (default = TRUE). This has no effect if not running
#'   interactively in RStudio.
#' @param force (logical) Force recreating of existing project database and
#'   config files (default = FALSE).
#' @param fetch_metadata (logical) Fetch each sample's NCBI metadata from its
#'   donor run and follow links to GBIF and GEOME (default = TRUE). Set FALSE
#'   when offline.
#' @param ... Additional arguments passed to [new_project_userAnnot()]
#'
#' @export
new_test_project_userAnnot <- function(
    path = ".",
    no_raw_data = FALSE,
    full_size = FALSE,
    executor = "local",
    container = paste0("macguigand/mitopilot:", utils::packageVersion("MitoPilot")),
    Rproj = TRUE,
    force = FALSE,
    fetch_metadata = TRUE,
    ...) {

  check_nextflow_version("new_test_project_userAnnot")

  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  path <- normalizePath(path)

  mapping <- app_sys(file.path("test_data", "mapping_test_userAnnot.csv")) |>
    utils::read.csv()

  dir.create(file.path(path, "assemblies"), showWarnings = FALSE)
  dir.create(file.path(path, "annotations"), showWarnings = FALSE)

  # Assemblies and annotations ----
  # Samples that reuse an assembly get their own copy under their ID.
  purrr::pwalk(mapping, function(...) {
    cur <- list(...)
    file.copy(app_sys(file.path("test_data", "assemblies", cur$Assembly)),
              file.path(path, "assemblies", paste0(cur$ID, ".fasta")), overwrite = TRUE)
    file.copy(app_sys(file.path("test_data", "annotations", cur$Annotation)),
              file.path(path, "annotations", cur$Annotation), overwrite = TRUE)
  })
  mapping$Assembly <- paste0(mapping$ID, ".fasta")

  # Reads ----
  if (!no_raw_data) {
    dir.create(file.path(path, "data"), showWarnings = FALSE)
    message("Fetching test data...")
    purrr::pwalk(mapping, function(...) {
      cur <- list(...)
      message(glue::glue("{cur$ID} - {cur$Taxon}"))
      purrr::walk(c("R1", "R2"), function(mate) {
        src <- fetch_test_reads(cur$Donors, mate, file.path(path, "data"), full_size)
        dest <- file.path(path, "data", glue::glue("{cur$ID}_{mate}.fastq.gz"))
        file.copy(src, dest, overwrite = TRUE)
      })
    })
  }

  mapping$BioSample <- mapping$Donors
  mapping$Donors <- NULL
  readr::write_csv(mapping, file.path(path, "mapping.csv"), quote = "none", na = "")

  new_project_userAnnot(
    path = path,
    mapping_fn = file.path(path, "mapping.csv"),
    mapping_id = "ID",
    data_path = if (!no_raw_data) file.path(path, "data", ""),
    no_raw_data = no_raw_data,
    assembly_path = file.path(path, "assemblies", ""),
    annotation_path = file.path(path, "annotations", ""),
    genetic_code = NULL,
    container = container,
    executor = executor,
    Rproj = Rproj,
    force = force,
    fetch_ncbi = fetch_metadata,
    link_sources = TRUE,
    ...
  )
}

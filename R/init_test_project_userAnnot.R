#' Initialize a test project for user-supplied annotations
#'
#' Sets up a seven-sample test project for the user-annotation workflow. Each
#' sample is one RefSeq teleost mitogenome from the GenBank records shipped with
#' MitoPilot: the record's sequence is the assembly (circular), its features
#' are the annotation, and small synthetic paired reads (about 20x) are sampled
#' from it. The coordinates are real; some names and extra features are
#' deliberately wrong so the gene-name review has work to do.
#'
#' \itemize{
#'   \item `GB_NCBI` (Danio rerio): an NCBI-style GFF3 with the full
#'     gene/mRNA/CDS/exon hierarchy. Every name is guessed.
#'   \item `GB_TBL` (Oncorhynchus mykiss): an NCBI 5-column feature table with
#'     anticodons. Every name is guessed.
#'   \item `GB_BED` (Gadus morhua): BED, with feature types inferred from the
#'     names and nad1 split into two blocks (a joined feature).
#'   \item `GB_TYPOS` (Thunnus thynnus): a GFF3 with odd names on the record's
#'     37 real genes. `mt-Nd2`, `MT-ND5`, `cytb`, `COII`, `COIII`,
#'     `tRNA-Leu (CUN)`, `trnS2(tga)`, and `l-rRNA` are guessed. `COX11` (cox1),
#'     `NDH4K` (nad4l), `ATP88` (atp8), a product-only
#'     `NADH dehydrogenase subnit 3` (nad3), `tRNA-Xxx` (trnK), `rrnX` (rrnS),
#'     `hypothetical protein` (nad6), and `ORF123` (atp6) are not. Each sits on
#'     a real gene, so mapping it by hand to the gene in brackets gives a valid
#'     annotation.
#'   \item `GB_DROP` (Xyrichtys novacula): a feature table with a `D-loop` and a
#'     control-region `misc_feature` (both guessed as the control region,
#'     `ctrl`), a light-strand `rep_origin` and a `repeat_region` (unmapped, so
#'     dropped), plus a `gap` and an `STS` (dropped on import), a
#'     pseudogene (nad6), ribosomal slippage (atp8 in two parts), a 5' partial
#'     nad1 with `codon_start` 2, a 3' partial rrnL, and nad5 spanning the
#'     origin. The two control-region features are adjacent, so validate
#'     reports `extra: ctrl` on purpose.
#'   \item `GB_MITOS` (Beryx splendens): a MITOS2-style GFF with lowercase
#'     names, numbered tRNAs (`trnL1(tag)`, `trnS2(tga)`), `OH` (guessed as
#'     `ctrl`), `OL` (unmapped), and cob spanning the origin as one row.
#'   \item `GB_GAPS` (Sebastes schlegelii): a plain GFF3 with correct names but
#'     five genes left out (cox2, nad3, trnD, trnK, rrnS). Turn on "Run
#'     annotation tools" in the Annotate Opts window and MitoPilot fills them in.
#' }
#'
#' @param path path to the directory for the test project (default = current
#'   working directory). Will be created if it does not already exist.
#' @param no_raw_data (logical) Set the project up without raw reads (default =
#'   FALSE). See [new_project_userAnnot()].
#' @param full_size (logical) Accepted for parity with [new_test_project_userAsmb()];
#'   has no effect, the reads are small synthetic files shipped with MitoPilot.
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
#' @param fetch_metadata (logical) Fetch NCBI metadata and follow links to GBIF
#'   and GEOME (default = TRUE). The demo samples have no BioSample, so nothing is
#'   fetched for them. Set FALSE when offline.
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
  file.copy(app_sys(file.path("test_data", "assemblies", mapping$Assembly)),
            file.path(path, "assemblies"), overwrite = TRUE)
  file.copy(app_sys(file.path("test_data", "annotations", mapping$Annotation)),
            file.path(path, "annotations"), overwrite = TRUE)

  if (!no_raw_data) {
    dir.create(file.path(path, "data"), showWarnings = FALSE)
    file.copy(app_sys(file.path("test_data", c(mapping$R1, mapping$R2))),
              file.path(path, "data"), overwrite = TRUE)
  }

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

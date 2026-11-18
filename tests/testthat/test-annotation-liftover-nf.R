asm_seq <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 3000, TRUE), collapse = ""))
set_at <- function(s, at, x) {
  substr(s, at, at + nchar(x) - 1L) <- x
  s
}
make_asm <- function(s = asm_seq, topology = "circular") {
  a <- Biostrings::DNAStringSet(s)
  names(a) <- paste("S1.1.1", topology)
  a
}
hit_row <- function(pos1, pos2, direction = "+") {
  data.frame(
    pos1 = as.integer(pos1), pos2 = as.integer(pos2), direction = direction,
    pctid = 97.14286, gap_pct = 1.428571, score = 60
  )
}
ref_row <- function(gene, type) {
  data.frame(gene = gene, type = type, pos1 = 1L, pos2 = 70L, direction = "+")
}
lift_args <- function(hit, ref, assembly = make_asm(), ref_dir = tempfile("noref")) {
  list(
    hit = hit, ref_feature = ref, accession = "NC_TEST.1",
    unit = list(ID = "S1", path = 1L, scaffold = 1L), assembly = assembly,
    coverage_fn = NULL, curate_params = params_fish_mito(), genetic_code = 2,
    ref_dir = ref_dir, template = data.frame(gene = character(), product = character())
  )
}
no_gap_hits <- function(...) {
  data.frame(
    acc = paste0("A", 1:3), Taxon = "t", eval = 0, target = "MKKK",
    pctid = 95, similarity = 95, gap_leading = 0L, gap_trailing = 0L
  )
}
stub_ref_dir <- function(gene) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(d, "featureProt"))
  writeLines(">ref\nMKKK", file.path(d, "featureProt", paste0(gene, ".fas")))
  d
}
# What liftover_curate.nf runs, in this process
nf_like_run <- function(st, ref_dir = "unused") {
  dir.create(file.path(st$dir, "curate"))
  dir.create(file.path(st$dir, "validate"))
  curate_mito_core(
    annotations_fn = st$ann_fn, assembly_fn = st$asm_fn, coverage_fn = NULL,
    genetic_code = st$genetic_code, out_dir = file.path(st$dir, "curate"),
    params = st$params_b64, ref_dir = ref_dir, feature_trim = FALSE,
    flip_rRNA_minus_strand = FALSE, ref_based_rc = FALSE
  )
  validate_mito_core(
    annotations_fn = file.path(st$dir, "curate", basename(st$ann_fn)),
    params = st$params_b64, out_dir = file.path(st$dir, "validate")
  )
  invisible(st$dir)
}

nf_state <- function(proj) {
  st <- do.call(liftover_prepare, c(
    lift_args(hit_row(1001, 1070), ref_row("rrnS", "rRNA")),
    dir = file.path(proj, ".liftover", "lift_test")
  ))
  st
}

test_that("only PCGs and rRNAs go through Nextflow", {
  expect_true(liftover_needs_nf("PCG"))
  expect_true(liftover_needs_nf("rRNA"))
  expect_false(liftover_needs_nf("tRNA"))
  expect_false(liftover_needs_nf(NA_character_))
})

test_that("prepare + finish reproduce liftover_annotate for a tRNA", {
  a <- lift_args(hit_row(1001, 1070), ref_row("trnF", "tRNA"))
  old <- do.call(liftover_annotate, a)
  st <- do.call(liftover_prepare, a)
  on.exit(unlink(st$dir, recursive = TRUE))
  expect_true(file.exists(st$ann_fn))
  expect_equal(read.csv(st$ann_fn)$gene, c("trnF", ".liftover_pad"))
  expect_null(liftover_curate_local(st))
  new <- liftover_finish(st, st$dir)
  expect_true(new$curated)
  drop_ts <- function(r) r[, setdiff(names(r), "time_stamp")]
  expect_equal(drop_ts(new$row), drop_ts(old$row))
})

test_that("finish with a skip reason gives the uncurated row and warning", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  st <- do.call(liftover_prepare, lift_args(hit_row(101, 400), ref_row("cox1", "PCG"), make_asm(s)))
  on.exit(unlink(st$dir, recursive = TRUE))
  res <- liftover_finish(st, NULL, skip_reason = "Docker is not running;\n")
  expect_false(res$curated)
  expect_match(res$row$warnings, "^Liftover: curation skipped \\(Docker is not running\\)")
  expect_no_match(res$row$warnings, "low reference similarity|check reference")
  expect_equal(c(res$row$start_codon, res$row$stop_codon), c("ATG", "TAA"))
  expect_equal(c(res$row$pos1, res$row$pos2), c(101L, 400L))
})

test_that("finish reads Nextflow-style outputs, keeping refHits and a lone T stop", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "T")
  st <- do.call(liftover_prepare, lift_args(hit_row(101, 398), ref_row("cox1", "PCG"), make_asm(s)))
  on.exit(unlink(st$dir, recursive = TRUE))
  testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = no_gap_hits,
    code = nf_like_run(st, stub_ref_dir("cox1"))
  )
  res <- liftover_finish(st, st$dir)
  r <- res$row
  expect_true(res$curated)
  expect_equal(nrow(r), 1L)
  expect_identical(r$stop_codon, "T")
  expect_equal(nchar(r$translation), 99L)
  expect_match(r$refHits, "A1")
  expect_equal(r$tool, "NC_TEST.1_liftover")
  expect_equal(r$edited, 1)
})

test_that("finish refuses outputs whose assembly changed", {
  st <- do.call(liftover_prepare, lift_args(hit_row(1001, 1070), ref_row("trnF", "tRNA")))
  on.exit(unlink(st$dir, recursive = TRUE))
  nf_like_run(st)
  a <- Biostrings::readDNAStringSet(st$asm_fn)
  Biostrings::writeXStringSet(Biostrings::reverseComplement(a), file.path(st$dir, "curate", basename(st$asm_fn)))
  expect_error(liftover_finish(st, st$dir), "Liftover aborted: curation changed the assembly")
})

test_that("a local reference directory is returned untouched", {
  proj <- withr::local_tempdir()
  ref <- withr::local_tempdir()
  dir.create(file.path(ref, "Clade", "featureProt"), recursive = TRUE)
  out <- liftover_ref_db_dir(ref, "Clade", proj)
  expect_equal(out, normalizePath(file.path(ref, "Clade")))
  expect_false(dir.exists(file.path(proj, ".liftover")))
  expect_equal(liftover_ref_db_dir(ref, "Clade.tar.gz", proj), out)
})

test_that("a local tarball is unpacked once into the project cache and reused", {
  proj <- withr::local_tempdir()
  ref <- withr::local_tempdir()
  src <- withr::local_tempdir()
  dir.create(file.path(src, "Clade", "featureProt"), recursive = TRUE)
  writeLines(">a\nMK", file.path(src, "Clade", "featureProt", "cox1.fas"))
  withr::with_dir(src, utils::tar(file.path(ref, "Clade.tar.gz"), "Clade", compression = "gzip"))
  out <- liftover_ref_db_dir(ref, "Clade.tar.gz", proj)
  expect_equal(out, normalizePath(file.path(proj, ".liftover", "ref_dbs", "Clade")))
  expect_true(file.exists(file.path(out, "featureProt", "cox1.fas")))
  writeLines("marker", file.path(out, "marker"))
  unlink(file.path(ref, "Clade.tar.gz"))
  expect_equal(liftover_ref_db_dir(ref, "Clade", proj), out)
  expect_true(file.exists(file.path(out, "marker")))
  expect_equal(list.files(file.path(proj, ".liftover", "ref_dbs")), "Clade")
})

test_that("a background fetch unpacks in a child process, then the cache is returned", {
  proj <- withr::local_tempdir()
  ref <- withr::local_tempdir()
  src <- withr::local_tempdir()
  dir.create(file.path(src, "Clade", "featureNuc"), recursive = TRUE)
  dir.create(file.path(src, "Clade", "featureProt"))
  writeLines(">a\nACGT", file.path(src, "Clade", "featureNuc", "rrnS.fas"))
  withr::with_dir(src, utils::tar(file.path(ref, "Clade.tar.gz"), "Clade", compression = "gzip"))
  p <- liftover_ref_db_dir(ref, "Clade", proj, background = TRUE)
  expect_s3_class(p, "process")
  expect_true(startsWith(basename(attr(p, "staging")), ".unpack_"))
  p$wait(60000)
  expect_equal(p$get_exit_status(), 0L)
  out <- liftover_ref_db_dir(ref, "Clade", proj, background = TRUE)
  expect_equal(out, normalizePath(file.path(proj, ".liftover", "ref_dbs", "Clade")))
  expect_true(file.exists(file.path(out, "featureNuc", "rrnS.fas")))
  expect_equal(list.files(file.path(proj, ".liftover", "ref_dbs"), all.files = TRUE, no.. = TRUE), "Clade")
  bad <- liftover_ref_db_dir("http://127.0.0.1:9/ref_dbs", "Other", proj, background = TRUE)
  bad$wait(60000)
  expect_equal(bad$get_exit_status(), 1L)
  expect_match(bad$read_all_output(), "^Could not download the curation reference database")
})

test_that("a tarball without the expected layout is refused and not cached", {
  proj <- withr::local_tempdir()
  ref <- withr::local_tempdir()
  src <- withr::local_tempdir()
  dir.create(file.path(src, "Other", "stuff"), recursive = TRUE)
  writeLines("x", file.path(src, "Other", "stuff", "a.txt"))
  withr::with_dir(src, utils::tar(file.path(ref, "Clade.tar.gz"), "Other", compression = "gzip"))
  expect_error(liftover_ref_db_dir(ref, "Clade", proj), "has no featureProt folder")
  expect_length(list.files(file.path(proj, ".liftover", "ref_dbs"), all.files = TRUE, no.. = TRUE), 0L)
  # a flat tarball (no top-level folder) with the right layout is accepted
  flat <- withr::local_tempdir()
  dir.create(file.path(flat, "featureProt"))
  writeLines(">a\nMK", file.path(flat, "featureProt", "cox1.fas"))
  withr::with_dir(flat, utils::tar(file.path(ref, "Flat.tar.gz"), "featureProt", compression = "gzip"))
  expect_true(dir.exists(file.path(liftover_ref_db_dir(ref, "Flat", proj), "featureProt")))
})

test_that("a fetch that loses the final rename to a complete cache succeeds", {
  d <- withr::local_tempdir()
  src <- withr::local_tempdir()
  dir.create(file.path(src, "Clade", "featureProt"), recursive = TRUE)
  withr::with_dir(src, utils::tar(file.path(d, "Clade.tar.gz"), "Clade", compression = "gzip"))
  out <- file.path(d, "cache", "Clade")
  dir.create(file.path(out, "featureProt"), recursive = TRUE)
  writeLines("x", file.path(out, "featureProt", "keep"))
  dir.create(file.path(d, "stage"))
  expect_equal(liftover_db_fetch("", file.path(d, "Clade.tar.gz"), file.path(d, "stage"), out, "Clade"), out)
  expect_true(file.exists(file.path(out, "featureProt", "keep")))
})

test_that("the fetch raises R's download timeout and restores it", {
  withr::local_options(timeout = 60)
  seen <- NULL
  testthat::with_mocked_bindings(
    .package = "utils",
    download.file = function(url, destfile, ...) {
      seen <<- getOption("timeout")
      stop("offline")
    },
    expect_error(liftover_db_fetch("http://x/y.tar.gz", tempfile(), tempfile(), tempfile(), "y"), "Could not download")
  )
  expect_equal(seen, 3600)
  expect_equal(getOption("timeout"), 60)
})

test_that("the image cache follows the pipeline unless the config or environment names one", {
  proj <- withr::local_tempdir()
  cfg <- file.path(proj, ".config")
  keys <- c("NXF_SINGULARITY_CACHEDIR", "NXF_APPTAINER_CACHEDIR")
  withr::local_envvar(NXF_SINGULARITY_CACHEDIR = NA, NXF_APPTAINER_CACHEDIR = NA)
  writeLines(c("singularity {", "  enabled = true", "}"), cfg)
  env <- liftover_sif_cache_env(cfg, proj)
  expect_equal(unname(env[keys]), rep(file.path(proj, "work", "singularity"), 2))
  withr::with_envvar(c(NXF_SINGULARITY_CACHEDIR = "/shared/sif"), expect_length(liftover_sif_cache_env(cfg, proj), 0L))
  expect_length(liftover_sif_cache_env(cfg, proj, c(NXF_APPTAINER_CACHEDIR = "/n/sif")), 0L)
  writeLines(c("singularity {", "  enabled = true", "  cacheDir = '/scratch/sif'", "}"), cfg)
  expect_length(liftover_sif_cache_env(cfg, proj), 0L)
  writeLines(c("apptainer.cacheDir = '/scratch/sif'"), cfg)
  expect_length(liftover_sif_cache_env(cfg, proj), 0L)
  writeLines(c("singularity {", "  // cacheDir = '/old'", "}"), cfg)
  expect_length(liftover_sif_cache_env(cfg, proj), 2L)
  st <- nf_state(proj)
  cmd <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    nf_pin_version = function(...) NA_character_,
    liftover_nf_command(st, file.path(proj, "db"), proj)
  )
  expect_equal(cmd$env[["NXF_SINGULARITY_CACHEDIR"]], file.path(proj, "work", "singularity"))
  expect_equal(cmd$env[["NXF_APPTAINER_CACHEDIR"]], file.path(proj, "work", "singularity"))
})

test_that("an unreachable database URL errors plainly and caches nothing", {
  proj <- withr::local_tempdir()
  expect_error(
    liftover_ref_db_dir("http://127.0.0.1:9/ref_dbs", "Clade", proj),
    "^Could not download the curation reference database"
  )
  expect_length(list.files(file.path(proj, ".liftover", "ref_dbs")), 0L)
  expect_error(liftover_ref_db_dir(withr::local_tempdir(), "Nope", proj), "not found")
})

test_that("the Nextflow command runs from the scratch folder under the project config", {
  proj <- withr::local_tempdir()
  cfg <- file.path(proj, ".config")
  writeLines(c("process {", "  executor = 'slurm'", "}"), cfg)
  st <- nf_state(proj)
  cmd <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    nf_pin_version = function(...) "25.10.6",
    liftover_nf_command(st, file.path(proj, "db"), proj)
  )
  a <- cmd$args
  expect_equal(cmd$command, "nextflow")
  expect_equal(cmd$wd, st$dir)
  expect_true(startsWith(st$dir, normalizePath(file.path(proj, ".liftover"))))
  expect_equal(a[which(a == "run") + 1], app_sys("nextflow", "liftover_curate.nf"))
  cs <- a[which(a == "-c") + 1]
  expect_equal(cs, c(cfg, file.path(st$dir, "liftover.config")))
  expect_equal(a[which(a == "-w") + 1], file.path(st$dir, "work"))
  expect_equal(a[which(a == "-log") + 1], file.path(st$dir, "nextflow.log"))
  expect_false("-resume" %in% a)
  expect_equal(a[which(a == "--lift_ref_db") + 1], file.path(proj, "db"))
  expect_equal(a[which(a == "--lift_outdir") + 1], st$dir)
  expect_equal(a[which(a == "--lift_params") + 1], st$params_b64)
  expect_equal(a[which(a == "--lift_coverage") + 1], app_sys("nextflow", "assets", "NO_FILE"))
  expect_equal(cmd$env[["NXF_VER"]], "25.10.6")
  expect_false("MITOPILOT_NO_CONDA" %in% names(cmd$env))
  expect_match(LIFTOVER_NF_CONFIG, "withName: 'liftover_curate'")
  expect_match(LIFTOVER_NF_CONFIG, "executor = 'local'")
  expect_match(LIFTOVER_NF_CONFIG, "memory = '2 GB'")
})

test_that("a native config puts its environment on the Nextflow launch", {
  skip_on_os("windows")
  proj <- withr::local_tempdir()
  act <- file.path(proj, "activate.sh")
  writeLines(c("export PATH=/opt/mp_native/bin:$PATH", "export MITOPILOT_NO_CONDA=1"), act)
  cfg <- file.path(proj, ".config")
  writeLines(paste0("params.native_activate = '", act, "'"), cfg)
  st <- nf_state(proj)
  cmd <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    native_nf_pin = function(nat_env) if (grepl("^/opt/mp_native/bin", nat_env[["PATH"]])) "25.10.4" else NA,
    liftover_nf_command(st, file.path(proj, "db"), proj)
  )
  expect_equal(cmd$env[["MITOPILOT_NO_CONDA"]], "1")
  expect_match(cmd$env[["PATH"]], "^/opt/mp_native/bin:")
  expect_equal(cmd$env[["NXF_VER"]], "25.10.4")
  expect_equal(cmd$env[[1]], "current")
})

test_that("a Nextflow failure is reduced to a plain reason", {
  log <- c(
    "N E X T F L O W  ~  version 25.10.6",
    "ERROR ~ Error executing process > 'liftover_curate'",
    "",
    "Caused by:",
    "  Process `liftover_curate` terminated with an error exit status (125)",
    "",
    "Command error:",
    "  Unable to find image 'nope/none:0' locally",
    "  docker: Error response from daemon: pull access denied for nope/none.",
    "",
    "Work dir:",
    "  /p/.liftover/lift_x/work/ab/cd"
  )
  expect_equal(
    liftover_nf_reason(log),
    "docker: Error response from daemon: pull access denied for nope/none."
  )
  expect_equal(liftover_nf_reason(c("ERROR ~ Unknown config attribute `foo`", "")), "Unknown config attribute `foo`")
  expect_equal(liftover_nf_reason(character()), "Nextflow stopped without a message")
  r <- c("Command error:", "  Error: Invalid JSON string.", "  Execution halted", "", "Work dir:")
  expect_equal(liftover_nf_reason(r), "Error: Invalid JSON string.")
})

ann_cols <- c(
  "ID", "path", "scaffold", "type", "gene", "product", "pos1", "pos2", "length",
  "direction", "anticodon", "tool", "start_codon", "stop_codon", "partial_start",
  "partial_stop", "translation", "notes", "warnings", "refHits", "edited", "time_stamp"
)

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
template0 <- data.frame(gene = character(), product = character())

lift <- function(hit, ref, assembly = make_asm(), ref_dir = tempfile("noref"),
                 template = template0, params = params_fish_mito(), coverage_fn = NULL) {
  liftover_annotate(
    hit = hit, ref_feature = ref, accession = "NC_TEST.1",
    unit = list(ID = "S1", path = 1L, scaffold = 1L), assembly = assembly,
    coverage_fn = coverage_fn, curate_params = params, genetic_code = 2,
    ref_dir = ref_dir, template = template
  )
}

stub_ref_dir <- function(gene) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  dir.create(file.path(d, "featureProt"))
  writeLines(">ref\nMKKK", file.path(d, "featureProt", paste0(gene, ".fas")))
  d
}
no_gap_hits <- function(...) {
  data.frame(
    acc = paste0("A", 1:3), Taxon = "t", eval = 0, target = "MKKK",
    pctid = 95, similarity = 95, gap_leading = 0L, gap_trailing = 0L
  )
}

test_that("tRNA liftover on a circular assembly keeps the hit and builds the row", {
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"))
  r <- res$row
  expect_equal(names(r), ann_cols)
  expect_equal(nrow(r), 1L)
  expect_equal(c(r$pos1, r$pos2, r$length), c(1001L, 1070L, 70L))
  expect_equal(r$direction, "+")
  expect_equal(r$tool, "NC_TEST.1_liftover")
  expect_true(is.na(r$anticodon))
  expect_equal(r$edited, 1)
  expect_equal(r$product, "tRNA-Phe")
  expect_equal(c(r$ID, r$path, r$scaffold), c("S1", "1", "1"))
  expect_match(r$notes, "^Liftover from NC_TEST\\.1: 97\\.1% identity, 1\\.4% gaps")
  expect_true(res$curated)
})

test_that("product comes from a template row of the same gene when present", {
  tmpl <- data.frame(gene = "trnF", product = "tRNA-Phenylalanine")
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), template = tmpl)
  expect_equal(res$row$product, "tRNA-Phenylalanine")
})

test_that("a minus-strand rRNA is not flipped", {
  res <- lift(hit_row(201, 1100, "-"), ref_row("rrnS", "rRNA"))
  expect_equal(res$row$direction, "-")
  expect_equal(c(res$row$pos1, res$row$pos2), c(201L, 1100L))
  expect_equal(res$row$product, "12S ribosomal RNA")
  expect_true(res$curated)
})

test_that("a hit in the middle of a linear assembly is not trimmed", {
  res <- lift(hit_row(1501, 1570), ref_row("trnF", "tRNA"), assembly = make_asm(topology = "linear"))
  expect_equal(c(res$row$pos1, res$row$pos2), c(1501L, 1570L))
  expect_true(res$curated)
})

test_that("a wrapping tRNA gets its true length", {
  p <- params_fish_mito()
  p$default_rules$tRNA$min_len <- 60
  res <- lift(hit_row(2971, 40), ref_row("trnF", "tRNA"), params = p)
  expect_equal(c(res$row$pos1, res$row$pos2), c(2971L, 40L))
  expect_equal(res$row$length, 70L)
  expect_false(grepl("below min length", res$row$warnings %|NA|% ""))
})

test_that("a PCG with an unreachable reference DB is added uncurated", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  res <- lift(hit_row(101, 400), ref_row("cox1", "PCG"), assembly = make_asm(s))
  r <- res$row
  expect_false(res$curated)
  expect_match(r$warnings, "Liftover: curation skipped (curation reference database not found)", fixed = TRUE)
  expect_no_match(r$warnings, "low reference similarity|check reference")
  expect_equal(c(r$start_codon, r$stop_codon), c("ATG", "TAA"))
  expect_equal(nchar(r$translation), 99L)
  expect_equal(r$product, "cytochrome c oxidase subunit 1")
})

test_that("a BLAST failure during curation degrades to an uncurated row", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  res <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = function(...) stop("blastp: command not found"),
    code = lift(hit_row(101, 400), ref_row("cox1", "PCG"),
                assembly = make_asm(s), ref_dir = stub_ref_dir("cox1"))
  )
  expect_false(res$curated)
  expect_match(res$row$warnings, "Liftover: curation skipped (blastp: command not found)", fixed = TRUE)
  expect_equal(c(res$row$pos1, res$row$pos2), c(101L, 400L))
})

test_that("a skip reason is kept on one line without semicolons", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  res <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = function(...) {
      warning("running command 'blastp' had status 127; see log")
      no_gap_hits()
    },
    code = lift(hit_row(101, 400), ref_row("cox1", "PCG"),
                assembly = make_asm(s), ref_dir = stub_ref_dir("cox1"))
  )
  expect_false(res$curated)
  expect_match(res$row$warnings, "Liftover: curation skipped (running command 'blastp' had status 127 see log)", fixed = TRUE)
})

test_that("an empty skip reason reads as unknown error", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  res <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = function(...) stop(""),
    code = lift(hit_row(101, 400), ref_row("cox1", "PCG"),
                assembly = make_asm(s), ref_dir = stub_ref_dir("cox1"))
  )
  expect_false(res$curated)
  expect_match(res$row$warnings, "Liftover: curation skipped (unknown error)", fixed = TRUE)
})

test_that("liftover runs with a per-base coverage file", {
  cov_fn <- withr::local_tempfile(fileext = ".csv")
  n <- nchar(asm_seq)
  write.csv(data.frame(
    SeqId = "S1.1.1", Position = seq_len(n), Call = strsplit(asm_seq, "")[[1]],
    Depth = 50L, Correct = 50L, ErrorRate = 0, MeanDepth = 50, GC = 0.5, mask = ""
  ), cov_fn, row.names = FALSE)
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), coverage_fn = cov_fn)
  expect_true(res$curated)
  expect_equal(c(res$row$pos1, res$row$pos2), c(1001L, 1070L))
  expect_false(grepl("low coverage|high error", res$row$warnings %|NA|% ""))
})

test_that("a curated PCG keeps a lone T stop codon as text", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "T")
  res <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = no_gap_hits,
    code = lift(hit_row(101, 398), ref_row("cox1", "PCG"),
                assembly = make_asm(s), ref_dir = stub_ref_dir("cox1"))
  )
  r <- res$row
  expect_true(res$curated)
  expect_identical(r$stop_codon, "T")
  expect_equal(nchar(r$translation), 99L)
  expect_false(grepl("non-standard stop codon", r$warnings %|NA|% ""))
  expect_false(is.na(r$refHits))
})

test_that("liftover aborts when curation changes the assembly", {
  expect_error(
    testthat::with_mocked_bindings(
      .package = "MitoPilot",
      curate_mito_core = function(annotations_fn, assembly_fn, out_dir, ...) {
        a <- Biostrings::readDNAStringSet(assembly_fn)
        Biostrings::writeXStringSet(Biostrings::reverseComplement(a), file.path(out_dir, basename(assembly_fn)))
        read.csv(annotations_fn)
      },
      code = lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"))
    ),
    "Liftover aborted: curation changed the assembly"
  )
})

test_that("liftover fails when curation drops the row", {
  expect_error(
    testthat::with_mocked_bindings(
      .package = "MitoPilot",
      curate_mito_core = function(annotations_fn, assembly_fn, out_dir, ...) {
        file.copy(assembly_fn, out_dir)
        read.csv(annotations_fn)[0, ]
      },
      code = lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"))
    ),
    "Liftover failed"
  )
})

asm_seq <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 3000, TRUE), collapse = ""))
set_at <- function(s, at, x) {
  substr(s, at, at + nchar(x) - 1L) <- x
  s
}
make_asm <- function(s = asm_seq) {
  a <- Biostrings::DNAStringSet(s)
  names(a) <- "S1.1.1 circular"
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
# An existing annotation as the app stores it
ex_row <- function(gene, type, pos1, pos2, direction = "+", warnings = NA_character_) {
  pcg <- type == "PCG"
  data.frame(
    ID = "S1", path = "1", scaffold = "1", type = type, gene = gene, product = NA_character_,
    pos1 = pos1, pos2 = pos2,
    length = if (pos1 > pos2) 3000 - pos1 + 1 + pos2 else pos2 - pos1 + 1,
    direction = direction, anticodon = NA_character_, tool = "mitos",
    start_codon = if (pcg) "ATG" else NA_character_, stop_codon = if (pcg) "TAA" else NA_character_,
    partial_start = 0, partial_stop = 0, translation = if (pcg) "MKKK" else NA_character_,
    notes = NA_character_, warnings = warnings, refHits = NA_character_, edited = 0,
    time_stamp = 1, fas = "nt", faa = if (pcg) "aa" else NA_character_
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
lift <- function(hit, ref, existing, assembly = make_asm(), params = params_fish_mito(),
                 ref_dir = tempfile("noref")) {
  liftover_annotate(
    hit = hit, ref_feature = ref, accession = "NC_TEST.1",
    unit = list(ID = "S1", path = 1L, scaffold = 1L), assembly = assembly,
    coverage_fn = NULL, curate_params = params, genetic_code = 2,
    ref_dir = ref_dir, template = data.frame(gene = character(), product = character()),
    existing = existing
  )
}
ws <- function(x) trimws(strsplit(x %|NA|% "", ";", fixed = TRUE)[[1]])

# trnF lifted to 1001-1070 overlaps cox1 (1066-1365) by 5 bp, under a quarter of either gene
test_that("a lifted tRNA gets overlap warnings from its neighbors", {
  p <- params_fish_mito()
  p$default_rules$tRNA$overlap <- list(start = 0, stop = FALSE)
  ex <- ex_row("cox1", "PCG", 1066, 1365)
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex, params = p)
  expect_true(res$curated)
  expect_true("exceeds max stop overlap" %in% ws(res$row$warnings))
  nb <- res$neighbors
  expect_equal(nrow(nb), 1L)
  expect_equal(c(nb$gene, nb$pos1), c("cox1", "1066"))
  expect_equal(ws(nb$warnings), "exceeds max start overlap")
})

test_that("a lifted gene with an existing name is a possible duplicate on both rows", {
  ex <- rbind(
    ex_row("trnF", "tRNA", 2501, 2570, warnings = "exceeds max length"),
    ex_row("trnP", "tRNA", 1801, 1870)
  )
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex)
  expect_true("possible duplicate" %in% ws(res$row$warnings))
  nb <- res$neighbors
  expect_equal(nrow(nb), 1L)
  expect_equal(c(nb$gene, nb$pos1), c("trnF", "2501"))
  expect_equal(nb$warnings, "exceeds max length; possible duplicate")
})

test_that("a warning an existing row already carries is not appended again", {
  ex <- ex_row("trnF", "tRNA", 2501, 2570, warnings = "possible duplicate")
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex)
  expect_true("possible duplicate" %in% ws(res$row$warnings))
  expect_equal(nrow(res$neighbors), 0L)
})

test_that("unaffected rows and warnings raised without the new gene are not neighbors", {
  ex <- rbind(
    ex_row("cox1", "PCG", 1066, 1365),
    ex_row("rrnS", "rRNA", 1501, 1900, warnings = "low coverage region"),
    ex_row("nad1", "PCG", 2001, 2300)
  )
  res <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex)
  nb <- res$neighbors
  expect_equal(nb$gene, "cox1")
  expect_equal(ws(nb$warnings), "exceeds max start overlap")
  expect_no_match(nb$warnings, "reference")
  expect_no_match(res$row$warnings %|NA|% "", "overlap")
})

test_that("no or zero-row existing annotations give the one-row result", {
  a <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), NULL)
  b <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex_row("trnF", "tRNA", 1, 70)[0, ])
  c <- lift(hit_row(1001, 1070), ref_row("trnF", "tRNA"), ex_row("trnF_DELETED_1", "tRNA", 2501, 2570))
  drop_ts <- function(r) r[, setdiff(names(r), "time_stamp")]
  expect_equal(drop_ts(b$row), drop_ts(a$row))
  expect_equal(drop_ts(c$row), drop_ts(a$row))
  expect_equal(nrow(a$neighbors), 0L)
  expect_equal(nrow(b$neighbors), 0L)
  expect_equal(nrow(c$neighbors), 0L)
  expect_equal(names(a$neighbors), c("gene", "pos1", "warnings"))
})

test_that("wrap-around existing and new rows validate with true lengths", {
  p <- params_fish_mito()
  p$default_rules$tRNA$min_len <- 60
  ex <- ex_row("trnP", "tRNA", 2961, 30)
  res <- lift(hit_row(2971, 40), ref_row("trnF", "tRNA"), ex, params = p)
  expect_equal(c(res$row$pos1, res$row$pos2, res$row$length), c(2971L, 40L, 70L))
  expect_false("below min length" %in% ws(res$row$warnings))
  expect_true("exceeds max overlap" %in% ws(res$row$warnings))
  nb <- res$neighbors
  expect_equal(nb$gene, "trnP")
  expect_equal(ws(nb$warnings), "exceeds max overlap")
})

test_that("an uncurated PCG still drops reference warnings but gets neighbor ones", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  ex <- ex_row("cox1", "PCG", 1501, 1800)
  res <- lift(hit_row(101, 400), ref_row("cox1", "PCG"), ex, assembly = make_asm(s))
  w <- ws(res$row$warnings)
  expect_false(res$curated)
  expect_true("Liftover: curation skipped (curation reference database not found)" %in% w)
  expect_true("possible duplicate" %in% w)
  expect_false(any(grepl("reference similarity|check reference", w)))
  expect_equal(res$neighbors$gene, "cox1")
  expect_equal(ws(res$neighbors$warnings), "possible duplicate")
})

test_that("a PCG curated through Nextflow outputs is validated against its neighbors", {
  s <- set_at(set_at(asm_seq, 101, "ATG"), 398, "TAA")
  st <- liftover_prepare(
    hit = hit_row(101, 400), ref_feature = ref_row("nad1", "PCG"), accession = "NC_TEST.1",
    unit = list(ID = "S1", path = 1L, scaffold = 1L), assembly = make_asm(s), coverage_fn = NULL,
    curate_params = params_fish_mito(), genetic_code = 2, ref_dir = "unused",
    template = data.frame(gene = character(), product = character())
  )
  on.exit(unlink(st$dir, recursive = TRUE))
  dir.create(file.path(st$dir, "curate"))
  testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = no_gap_hits,
    code = curate_mito_core(
      annotations_fn = st$ann_fn, assembly_fn = st$asm_fn, coverage_fn = NULL,
      genetic_code = 2, out_dir = file.path(st$dir, "curate"), params = st$params_b64,
      ref_dir = stub_ref_dir("nad1"), feature_trim = FALSE,
      flip_rRNA_minus_strand = FALSE, ref_based_rc = FALSE
    )
  )
  ex <- rbind(ex_row("nad1", "PCG", 1501, 1800), ex_row("trnF", "tRNA", 391, 460))
  res <- liftover_finish(st, st$dir, existing = ex)
  w <- ws(res$row$warnings)
  expect_true(res$curated)
  expect_true("possible duplicate" %in% w)
  expect_true("exceeds max stop overlap" %in% w)
  expect_false(any(grepl("reference", w)))
  expect_equal(res$neighbors$gene, "nad1")
})

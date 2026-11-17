test_that("curate_db_is_custom separates the bundled database from anything else", {
  expect_false(curate_db_is_custom("Metazoa_RefSeq235"))
  expect_true(curate_db_is_custom("Metazoa_RefSeq235_custom"))
  expect_false(curate_db_is_custom(NULL))
  expect_false(curate_db_is_custom(NA_character_))
})

test_that("curate_db_fields resolves the form to ref_dir / ref_db", {
  expect_equal(
    curate_db_fields(FALSE, "ignored", "https://x/Mitos2"),
    list(ref_dir = "https://x/Mitos2", ref_db = "Metazoa_RefSeq235")
  )

  td <- tempfile()
  db <- file.path(td, "Metazoa_RefSeq235_custom")
  dir.create(file.path(db, "featureProt"), recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE))
  ref <- curate_db_fields(TRUE, paste0(" ", db, "/ "), "unused")
  expect_equal(ref, list(ref_dir = td, ref_db = "Metazoa_RefSeq235_custom"))
  # round trip: the stored pair rebuilds the path shown in the form
  expect_equal(file.path(ref$ref_dir, ref$ref_db), db)
  expect_true(curate_db_is_custom(ref$ref_db))

  expect_match(curate_db_fields(TRUE, td, "unused")$error, "featureProt")
  expect_match(curate_db_fields(TRUE, "", "unused")$error, "featureProt")
})

test_that("custom_curation_db only builds on the Metazoa base", {
  expect_error(custom_curation_db(base_db = "Chordata"), "Metazoa_RefSeq235")
})

test_that("get_top_hits errors when BLAST fails instead of returning no hits", {
  skip_if(Sys.which("blastp") == "", "blastp not installed")
  fa <- tempfile(fileext = ".fas")
  writeLines(c(">a test", "MAHPLQLGFQDAMSPIMEELLHFHDHTLMIVF"), fa)
  on.exit(unlink(fa))
  # no BLAST index beside the FASTA, so blastp exits non-zero
  expect_error(
    suppressWarnings(get_top_hits(fa, "MAHPLQLGFQDAMSPIMEEL", condaenv = NULL)),
    "blastp failed"
  )
})

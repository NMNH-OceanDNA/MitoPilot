# validate_mito_core() treats rows sharing a JOIN group as one gene copy.

jL <- 2000L

jvrow <- function(pos1, pos2, notes = NA_character_, start_codon = "ATG",
                  stop_codon = "TAA", translation = strrep("M", 60L),
                  join_translation = NULL) {
  out <- data.frame(
    contig = "ctg1", type = "PCG", gene = "cox1", product = "cox1",
    pos1 = as.integer(pos1), pos2 = as.integer(pos2),
    length = as.integer(pos2 - pos1 + 1L), direction = "+",
    start_codon = start_codon, stop_codon = stop_codon, anticodon = NA_character_,
    translation = translation, refHits = "{}", notes = notes,
    stringsAsFactors = FALSE
  )
  if (!is.null(join_translation)) out$join_translation <- join_translation
  out
}

run_jvalidate <- function(ann) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  ann_fn <- file.path(d, "x_annotations_1.csv")
  write.csv(ann, ann_fn, row.names = FALSE)
  cov_fn <- file.path(d, "cov.csv")
  write.csv(
    data.frame(
      SeqId = "ctg1", Position = seq_len(jL), Call = "A",
      MeanDepth = 100, GC = 0.4, ErrorRate = 0
    ),
    cov_fn, row.names = FALSE
  )
  validate_ctenophore_mito(
    annotations_fn = ann_fn, coverage_fn = cov_fn,
    params = params_ctenophore_mito(), out_dir = d
  )
}

test_that("a JOIN group counts once toward the gene count", {
  res <- run_jvalidate(rbind(
    jvrow(100, 400, "JOIN: mode=exon group=1"),
    jvrow(600, 900, "JOIN: mode=exon group=1")
  ))
  expect_false(any(grepl("possible duplicate", res$annotations$warnings)))
  expect_true(is.na(res$summary$extra) || !grepl("cox1", res$summary$extra))
  expect_false(isTRUE(grepl("cox1", res$summary$missing)))
  expect_equal(res$summary$PCGCount, 1L)
  expect_equal(lengths(regmatches(res$summary$structure,
                                  gregexpr("cox1", res$summary$structure))), 1L)
})

test_that("the same rows without a JOIN note are still extra", {
  res <- run_jvalidate(rbind(jvrow(100, 400), jvrow(600, 900)))
  expect_match(res$summary$extra, "cox1")
  expect_true(all(grepl("possible duplicate", res$annotations$warnings)))
})

test_that("codon checks read the outer members and the spliced translation", {
  res <- run_jvalidate(rbind(
    jvrow(600, 900, "JOIN: mode=exon group=1", start_codon = "CCC",
          stop_codon = "TAA", translation = "MK*F", join_translation = "MKFL"),
    jvrow(100, 400, "JOIN: mode=exon group=1", start_codon = "ATG",
          stop_codon = "CCC", translation = "MK*F", join_translation = "MKFL")
  ))
  w <- dplyr::coalesce(res$annotations$warnings, "")
  expect_false(any(grepl("non-standard start codon", w)))
  expect_false(any(grepl("non-standard stop codon", w)))
  expect_false(any(grepl("internal stop codon", w)))
  expect_false("join_translation" %in% names(res$annotations))

  res <- run_jvalidate(rbind(
    jvrow(100, 400, "JOIN: mode=exon group=1", start_codon = "CCC",
          translation = "MKF", join_translation = "MK*F"),
    jvrow(600, 900, "JOIN: mode=exon group=1", stop_codon = "CCC",
          translation = "MKF", join_translation = "MK*F")
  ))
  w <- dplyr::coalesce(res$annotations$warnings, "")
  expect_match(w[1], "non-standard start codon")
  expect_false(grepl("non-standard start codon", w[2]))
  expect_match(w[2], "non-standard stop codon")
  expect_false(grepl("non-standard stop codon", w[1]))
  expect_equal(sum(grepl("internal stop codon", w)), 1L)
})

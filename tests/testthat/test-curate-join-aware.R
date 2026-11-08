# Pure parts of join-aware curation (the curate loop itself needs BLAST).

join_seq <- function() Biostrings::DNAString("GGATGAACCCATTTTAAGGG")

jrow <- function(pos1, pos2, direction) {
  data.frame(
    pos1 = pos1, pos2 = pos2, direction = direction,
    start_codon = "ATG", stop_codon = "TAA",
    partial_start = 0L, partial_stop = 0L,
    stringsAsFactors = FALSE
  )
}

test_that("is_join_member flags rows with a leading JOIN marker", {
  expect_equal(
    is_join_member(c("JOIN: mode=exon group=1", "JOIN: mode=frameshift group=2 note=x; y",
                     NA, "", "x; JOIN: mode=exon group=1")),
    c(TRUE, TRUE, FALSE, FALSE, FALSE)
  )
})

test_that("splice_join_cds splices a plus-strand group", {
  gc <- Biostrings::getGeneticCode("2")
  m <- rbind(jrow(11, 17, "+"), jrow(3, 7, "+"))
  expect_equal(splice_join_cds(m, join_seq(), gc)$translation, "MKF")
})

test_that("splice_join_cds splices a minus-strand group", {
  gc <- Biostrings::getGeneticCode("2")
  rc <- Biostrings::reverseComplement(join_seq())
  m <- rbind(jrow(4, 10, "-"), jrow(14, 18, "-"))
  expect_equal(splice_join_cds(m, rc, gc)$translation, "MKF")
})

test_that("curate BLASTs a JOIN group once and leaves its members in place", {
  s <- paste(rep("AAC", 300L), collapse = "")
  substr(s, 100, 102) <- "ATG"
  substr(s, 598, 600) <- "TAA"
  d <- withr::local_tempdir()
  dir.create(file.path(d, "featureProt"))
  writeLines(">ref\nMKKK", file.path(d, "featureProt", "nad1.fas"))
  asm <- Biostrings::DNAStringSet(s)
  names(asm) <- "ctg1 topology=linear"
  Biostrings::writeXStringSet(asm, file.path(d, "assembly.fasta"))
  ann <- data.frame(
    contig = "ctg1", type = c("PCG", "PCG", "tRNA"), gene = c("nad1", "nad1", "trnF"),
    product = "x", pos1 = c(100L, 400L, 600L), pos2 = c(300L, 600L, 670L),
    length = 0L, direction = "+", start_codon = c("ATG", "AAC", NA),
    stop_codon = c("AAC", "TAA", NA), translation = c("MNN", "NNN", NA),
    notes = c("JOIN: mode=exon group=1", "JOIN: mode=exon group=1", NA),
    stringsAsFactors = FALSE
  )
  write.csv(ann, file.path(d, "annotations.csv"), row.names = FALSE)
  dir.create(file.path(d, "out"))

  seen <- character(0)
  res <- testthat::with_mocked_bindings(
    .package = "MitoPilot",
    get_top_hits = function(db, translation, ...) {
      seen <<- c(seen, translation)
      data.frame(acc = "A1", Taxon = "t", eval = 0, target = "nad1", pctid = 90,
                 similarity = 90, gap_leading = 2L, gap_trailing = 0L)
    },
    code = curate_ctenophore_mito(
      annotations_fn = file.path(d, "annotations.csv"),
      assembly_fn = file.path(d, "assembly.fasta"),
      out_dir = file.path(d, "out"),
      params = params_ctenophore_mito(),
      ref_dir = d, genetic_code = 4, feature_trim = FALSE
    )
  ) |> as.data.frame()

  expect_equal(seen, paste0("M", strrep("N", 132L)))
  pcg <- res[res$type == "PCG", ]
  expect_equal(pcg$pos1, c(100L, 400L))
  expect_equal(pcg$pos2, c(300L, 600L))
  expect_equal(pcg$stop_codon[2], "TAA")
  expect_equal(pcg$refHits[1], pcg$refHits[2])
  expect_equal(pcg$join_translation, rep(seen, 2L))
})

test_that("remote reference prepend gives every JOIN member the same refHits", {
  d <- withr::local_tempdir()
  seq <- paste0("ATG", strrep("AAC", 132), "TAA")
  f <- file.path(d, "ref.json")
  jsonlite::write_json(list(
    accession = "X1", organism = "o", sequence = seq, genetic_code = 4,
    pcg = data.frame(gene = "nad1", pos1 = 1, pos2 = nchar(seq), direction = "+")
  ), f, auto_unbox = TRUE)
  shared <- json_string(data.frame(
    acc = "A1", Taxon = "t", eval = 0, target = "x", pctid = 90,
    similarity = 90, gap_leading = 0L, gap_trailing = 0L
  ))
  ann <- data.frame(
    type = "PCG", gene = "nad1", translation = c("MNN", "NNN"),
    join_translation = paste0("M", strrep("N", 132)), refHits = shared,
    stringsAsFactors = FALSE
  )
  out <- suppressMessages(prepend_blast_hit_to_refhits(ann, f))
  expect_identical(out$refHits[1], out$refHits[2])
  top <- json_parse(out$refHits[1], TRUE)
  expect_equal(top$acc[1], "X1")
  expect_equal(top$gap_leading[1], 0)
})

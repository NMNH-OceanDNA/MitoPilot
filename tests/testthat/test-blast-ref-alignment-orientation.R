test_that("a forward alignment passes through untouched", {
  aln <- data.frame(
    aligned_sample = "ACGTACGT", aligned_ref = "ACGTACGT",
    rotation = 0L, ref_length = 8L, ref_start = 0L, strand = "+",
    stringsAsFactors = FALSE
  )
  out <- normalize_blast_ref_alignment(aln)
  expect_equal(out$aligned_sample, aln$aligned_sample)
  expect_equal(out$aligned_ref, aln$aligned_ref)
  expect_equal(out$strand, "+")
  expect_false(out$ref_rc)
})

test_that("a reverse alignment is flipped into the sample's orientation", {
  # compute_blast_ref_alignment() stores the sample reverse-complemented when that
  # strand scored higher, which mirrors the synteny plot against the coverage map.
  # Normalising puts the sample back in its stored frame and moves the flip to the
  # reference.
  sample_stored <- "AAACCCGGGTTT"
  ref           <- "AAACCCGGGTTA"
  aln <- data.frame(
    aligned_sample = as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(sample_stored))),
    aligned_ref = ref,
    rotation = 0L, ref_length = nchar(ref), ref_start = 0L, strand = "-",
    stringsAsFactors = FALSE
  )
  out <- normalize_blast_ref_alignment(aln)
  expect_equal(out$aligned_sample, sample_stored)
  expect_equal(out$aligned_ref,
               as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref))))
  expect_equal(out$strand, "+")
  expect_true(out$ref_rc)
})

test_that("gap columns survive the flip and stay paired", {
  s <- "AC-GTTGCA"
  r <- "ACGGT-GCA"
  aln <- data.frame(
    aligned_sample = s, aligned_ref = r, rotation = 0L,
    ref_length = 8L, ref_start = 0L, strand = "-", stringsAsFactors = FALSE
  )
  out <- normalize_blast_ref_alignment(aln)
  # Same number of columns, so every position still pairs with its partner.
  expect_equal(nchar(out$aligned_sample), nchar(s))
  expect_equal(nchar(out$aligned_ref), nchar(r))
  # Gaps are complemented to themselves, and land at the mirrored column.
  expect_equal(gregexpr("-", out$aligned_sample)[[1]][1],
               nchar(s) - gregexpr("-", s)[[1]][1] + 1L)
  expect_equal(gregexpr("-", out$aligned_ref)[[1]][1],
               nchar(r) - gregexpr("-", r)[[1]][1] + 1L)
})

test_that("normalisation tolerates empty and malformed input", {
  expect_null(normalize_blast_ref_alignment(NULL))
  empty <- data.frame(aligned_sample = character(), aligned_ref = character(),
                      strand = character(), stringsAsFactors = FALSE)
  expect_equal(nrow(normalize_blast_ref_alignment(empty)), 0L)
  # No strand column (pre-v2 aligner output): flagged not-flipped, left alone.
  old <- data.frame(aligned_sample = "ACGT", aligned_ref = "ACGT",
                    stringsAsFactors = FALSE)
  out <- normalize_blast_ref_alignment(old)
  expect_false(out$ref_rc)
  expect_equal(out$aligned_sample, "ACGT")
})

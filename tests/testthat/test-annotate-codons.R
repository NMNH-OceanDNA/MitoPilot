test_that("add_pcg_codons fills codons for PCG rows and matches header first token", {
  asm <- Biostrings::DNAStringSet(c("S.1.1 circular" = paste0("ATG", "AAACCC", "TAA", "GGGG")))
  ann <- tibble::tibble(contig = "S.1.1", type = c("PCG", "tRNA"), gene = c("cox1", "trnF"),
                        pos1 = c(1L, 13L), pos2 = c(12L, 16L), direction = "+")
  out <- add_pcg_codons(ann, asm, genetic_code = "2")
  expect_equal(out$start_codon, c("ATG", NA))
  expect_equal(out$stop_codon, c("TAA", NA))
  expect_equal(out$translation, c("MKP", NA))
  expect_equal(out$gene, ann$gene)
})

test_that("add_pcg_codons handles minus strand and origin wrap", {
  seq <- paste0("GGGTTTCAT", "ACACACAC", "TTA")
  asm <- Biostrings::DNAStringSet(c("S.1.1 circular" = seq))
  ann <- tibble::tibble(contig = "S.1.1", type = "PCG", gene = "cox1",
                        pos1 = 18L, pos2 = 9L, direction = "-")
  out <- add_pcg_codons(ann, asm, genetic_code = "2")
  cds <- Biostrings::reverseComplement(Biostrings::DNAString(paste0("TTA", "GGGTTTCAT")))
  expect_equal(out$start_codon, "ATG")
  expect_equal(out$stop_codon, "TAA")
  expect_equal(
    out$translation,
    as.character(Biostrings::translate(cds[1:9], genetic.code = Biostrings::getGeneticCode("2")))
  )
})

test_that("add_pcg_codons honors codon_start on 5'-partial CDS", {
  mk <- function(seq, len) {
    asm <- Biostrings::DNAStringSet(c("S.1.1" = seq))
    ann <- tibble::tibble(contig = "S.1.1", type = "PCG", gene = "cox1",
                          pos1 = 1L, pos2 = len, direction = "+", codon_start = 2L)
    add_pcg_codons(ann, asm, genetic_code = "2")
  }
  full <- mk("GATGAAACCCTAAGG", 13L)
  expect_equal(full$stop_codon, "TAA")
  expect_equal(full$translation, "MKP")
  expect_equal(full$start_codon, "GAT")
  part <- mk("GATGAAACCCTGG", 11L)
  expect_equal(part$stop_codon, "T")
  expect_equal(part$translation, "MKP")
  asm <- Biostrings::DNAStringSet(c("S.1.1" = "ATGAAACCCTAAGGGG"))
  ann <- tibble::tibble(contig = "S.1.1", type = "PCG", gene = "cox1",
                        pos1 = 1L, pos2 = 12L, direction = "+", codon_start = NA_integer_)
  expect_equal(add_pcg_codons(ann, asm, "2")$translation, "MKP")
})

# Curation can extend the start gene's start back across the origin of a
# circular contig; the contig is then re-rotated so the gene starts at 1 again.

reanchor_fixture <- function(name = "ctg1 circular") {
  seq_chr <- paste(c("A", "C", "G", "T")[(seq_len(40) * 7L) %% 4L + 1L], collapse = "")
  asm <- Biostrings::DNAStringSet(seq_chr)
  names(asm) <- name
  ann <- data.frame(
    contig = "ctg1", gene = c("nad1", "cox1"),
    pos1 = c(15L, 35L), pos2 = c(20L, 10L), direction = "+"
  )
  cov <- data.frame(SeqId = "ctg1", Position = 1:40, Call = strsplit(seq_chr, "")[[1]])
  list(asm = asm, ann = ann, cov = cov, seq = seq_chr, key = c(ctg1 = name))
}

test_that("an origin-spanning start gene is rotated back to position 1", {
  f <- reanchor_fixture()
  expect_message(
    res <- .reanchor_start_gene(f$asm, f$ann, f$cov, f$key, c(ctg1 = "cox1")),
    "re-rotated"
  )
  expect_equal(as.character(res$assembly[[1]]),
               paste0(substr(f$seq, 35, 40), substr(f$seq, 1, 34)))
  expect_equal(names(res$assembly), "ctg1 circular")
  expect_equal(res$annotations$gene, c("cox1", "nad1"))
  expect_equal(res$annotations$pos1, c(1L, 21L))
  expect_equal(res$annotations$pos2, c(16L, 26L))
  expect_equal(res$coverage$Position, 1:40)
  expect_equal(paste(res$coverage$Call, collapse = ""), as.character(res$assembly[[1]]))
})

test_that("contigs are left alone when there is nothing to fix", {
  f <- reanchor_fixture()
  same <- function(res) {
    expect_equal(as.character(res$assembly[[1]]), f$seq)
    expect_equal(res$annotations, f$ann)
  }
  # linear contig
  lin <- reanchor_fixture("ctg1 linear")
  same(.reanchor_start_gene(lin$asm, lin$ann, lin$cov, lin$key, c(ctg1 = "cox1")))
  # anchor gone, or flipped onto the minus strand
  same(.reanchor_start_gene(f$asm, f$ann, f$cov, f$key, c(ctg1 = "atp6")))
  f$ann$direction <- "-"
  same(.reanchor_start_gene(f$asm, f$ann, f$cov, f$key, c(ctg1 = "cox1")))
  # no anchor recorded
  same(.reanchor_start_gene(f$asm, f$ann, NULL, f$key, character(0)))
})

ua_seq <- paste0(
  strrep("ACGT", 7), "AC",                   # 1-30
  strrep("GGCA", 17), "GG",                  # 31-100 trnF
  "ATG", strrep("GCT", 98), "TAA",           # 101-400 cox1
  strrep("ACGT", 50)                         # 401-600
)

ua_fixture <- function(topology = "linear", pos2_cox1 = 400L) {
  d <- tempfile("ua_in")
  dir.create(d)
  fa <- file.path(d, "S_assembly_1.fasta")
  asm <- Biostrings::DNAStringSet(ua_seq)
  names(asm) <- paste("S.1.1", topology)
  Biostrings::writeXStringSet(asm, fa)
  n <- nchar(ua_seq)
  cov <- file.path(d, "S_coverageStats_1.csv")
  utils::write.csv(data.frame(
    SeqId = "S.1.1", Position = seq_len(n), Call = strsplit(ua_seq, "")[[1]],
    Depth = 40, Correct = 40, ErrorRate = 0, MeanDepth = 40, GC = 0.5
  ), cov, row.names = FALSE)
  csv <- file.path(d, "S_user_annotations.csv")
  utils::write.csv(data.frame(
    contig = "S.1.1", type = c("tRNA", "PCG"), gene = c("trnF", "cox1"),
    product = c("tRNA-Phe", "cytochrome c oxidase subunit I"),
    pos1 = c(31L, 101L), pos2 = c(100L, pos2_cox1), length = c(70L, 300L),
    direction = "+", anticodon = c("gaa", NA), partial_start = c(NA, 1L),
    partial_stop = c(NA, 0L), codon_start = c(NA, 1L), notes = c(NA, "imported"),
    tool = "user"
  ), csv, row.names = FALSE, na = "")
  out <- tempfile("ua_out")
  dir.create(out)
  list(fa = fa, cov = cov, csv = csv, out = out)
}

ua_run <- function(f, ...) {
  annotate(
    assembly_fn = f$fa, coverage_fn = f$cov, cpus = 1, genetic_code = "2",
    ref_db = "x", ref_dir = ".", use_mitos = FALSE, use_trnaScan = FALSE,
    use_arwen = FALSE, use_aragorn = FALSE, use_mitofinder = FALSE,
    out_dir = f$out, ...
  )
}

test_that("annotate() with user_annotations_fn skips tools and emits codons", {
  f <- ua_fixture()
  ua_run(f, start_gene = NA, coverage_trim = FALSE, user_annotations_fn = f$csv)
  res <- read.csv(list.files(f$out, "_annotations_", full.names = TRUE))
  expect_equal(sort(res$gene), c("cox1", "trnF"))
  expect_false(is.na(res$translation[res$gene == "cox1"]))
  expect_equal(res$start_codon[res$gene == "cox1"], "ATG")
  expect_equal(res$stop_codon[res$gene == "cox1"], "TAA")
  expect_equal(unique(res$tool), "user")
  expect_true(all(c("partial_start", "partial_stop", "codon_start", "notes") %in% names(res)))
  expect_equal(res$partial_start[res$gene == "cox1"], 1L)
  expect_equal(res$notes[res$gene == "cox1"], "imported")
})

test_that("annotate() errors when no annotation source is active", {
  f <- ua_fixture()
  expect_error(ua_run(f, start_gene = NA, coverage_trim = FALSE), "no annotation source")
})

test_that("annotate() errors when an imported feature exceeds the contig", {
  f <- ua_fixture(pos2_cox1 = nchar(ua_seq) + 10L)
  expect_error(
    ua_run(f, start_gene = NA, coverage_trim = FALSE, user_annotations_fn = f$csv),
    "1 feature\\(s\\) extend past the end of contig S.1.1"
  )
})

test_that("start_gene = NA skips rotation without warning", {
  f <- ua_fixture("circular")
  expect_no_warning(
    ua_run(f, start_gene = NA, coverage_trim = FALSE, user_annotations_fn = f$csv)
  )
  res <- read.csv(list.files(f$out, "_annotations_", full.names = TRUE))
  expect_equal(res$pos1[res$gene == "trnF"], 31L)
})

test_that("coverage_trim is ignored for imported annotations", {
  f <- ua_fixture()
  expect_message(
    ua_run(f, start_gene = NA, coverage_trim = TRUE, user_annotations_fn = f$csv),
    "coverage_trim ignored"
  )
})

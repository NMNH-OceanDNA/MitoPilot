# assess_partial_ends(): curation owns partial_start / partial_stop.

partial_asm <- function() {
  asm <- Biostrings::DNAStringSet(c(strrep("A", 1000L), strrep("A", 1000L)))
  names(asm) <- c("L.1.1 linear", "C.1.1 circular")
  asm
}

partial_key <- function(asm) setNames(names(asm), stringr::str_extract(names(asm), "^\\S+"))

partial_rules <- function() {
  cod <- list(start_codons = c("ATG", "GTG"), stop_codons = c("TAA", "TAG", "T", "TA"))
  list(cox1 = cod, nad5 = cod, nad6 = cod, rrnL = cod, trnF = cod)
}

prow <- function(gene, type, pos1, pos2, direction = "+", contig = "L.1.1",
                 start_codon = NA_character_, stop_codon = NA_character_,
                 partial_start = 0L, partial_stop = 0L, notes = NA_character_) {
  tibble::tibble(
    contig = contig, type = type, gene = gene,
    pos1 = as.integer(pos1), pos2 = as.integer(pos2), direction = direction,
    start_codon = start_codon, stop_codon = stop_codon,
    partial_start = partial_start, partial_stop = partial_stop, notes = notes
  )
}

run_partial <- function(ann) {
  asm <- partial_asm()
  assess_partial_ends(ann, asm, partial_key(asm), partial_rules())
}

test_that("edge of a linear contig marks the strand-correct end partial", {
  out <- run_partial(dplyr::bind_rows(
    prow("trnF", "tRNA", 1, 70),
    prow("rrnL", "rRNA", 100, 1000, direction = "-"),
    prow("nad6", "PCG", 1, 300, direction = "-", start_codon = "ATG", stop_codon = "TAA")
  ))
  expect_equal(out$partial_start[out$gene == "trnF"], 1L)
  expect_equal(out$partial_stop[out$gene == "trnF"], 0L)
  expect_equal(out$partial_start[out$gene == "rrnL"], 1L)
  expect_equal(out$partial_stop[out$gene == "rrnL"], 0L)
  expect_equal(out$partial_start[out$gene == "nad6"], 0L)
  expect_equal(out$partial_stop[out$gene == "nad6"], 1L)
})

test_that("PCG codon evidence: bad start -> 5' partial, bad full stop -> 3' partial, short stop -> complete", {
  out <- run_partial(dplyr::bind_rows(
    prow("cox1", "PCG", 200, 400, start_codon = "CTG", stop_codon = "TAA"),
    prow("cox1", "PCG", 401, 600, start_codon = "ATG", stop_codon = "AAA"),
    prow("cox1", "PCG", 601, 700, start_codon = "ATG", stop_codon = "TA"),
    prow("orfX", "PCG", 701, 800, start_codon = "CTG", stop_codon = "AAA")
  ))
  expect_equal(out$partial_start[out$gene == "cox1" & out$start_codon == "CTG"], 1L)
  expect_equal(out$partial_stop[out$gene == "cox1" & out$start_codon == "CTG"], 0L)
  expect_equal(out$partial_stop[out$gene == "cox1" & out$stop_codon == "AAA"], 1L)
  expect_equal(out$partial_start[out$gene == "cox1" & out$stop_codon == "AAA"], 0L)
  expect_equal(out$partial_stop[out$gene == "cox1" & out$stop_codon == "TA"], 0L)
  # no rule for the gene: no codon evidence
  expect_equal(out$partial_start[out$gene == "orfX"], 0L)
  expect_equal(out$partial_stop[out$gene == "orfX"], 0L)
})

test_that("circular contig has no edges", {
  out <- run_partial(dplyr::bind_rows(
    prow("cox1", "PCG", 1, 1000, contig = "C.1.1", start_codon = "ATG", stop_codon = "TAA"),
    prow("trnF", "tRNA", 1, 70, direction = "-", contig = "C.1.1")
  ))
  expect_equal(out$partial_start[out$contig == "C.1.1"], c(0L, 0L))
  expect_equal(out$partial_stop[out$contig == "C.1.1"], c(0L, 0L))
})

test_that("JOIN group: only outer ends assessed", {
  out <- run_partial(dplyr::bind_rows(
    prow("nad5", "PCG", 600, 800, start_codon = "ATG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1"),
    prow("nad5", "PCG", 100, 300, start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1")
  ))
  i_5p <- which(out$pos1 == 100)
  i_3p <- which(out$pos1 == 600)
  expect_equal(out$partial_start[i_5p], 1L)
  expect_equal(out$partial_stop[i_5p], 0L)
  expect_equal(out$partial_start[i_3p], 0L)
  expect_equal(out$partial_stop[i_3p], 1L)

  # minus strand: 5' member has the largest pos2; group crossing a circular origin
  out <- run_partial(dplyr::bind_rows(
    prow("nad5", "PCG", 100, 300, direction = "-", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1"),
    prow("nad5", "PCG", 600, 800, direction = "-", start_codon = "ATG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1"),
    prow("cox1", "PCG", 900, 50, contig = "C.1.1", start_codon = "ATG", stop_codon = "TAA",
         notes = "JOIN: mode=exon group=2"),
    prow("cox1", "PCG", 80, 200, contig = "C.1.1", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=2")
  ))
  expect_equal(out$partial_start, c(0L, 0L, 0L, 0L))
  expect_equal(out$partial_stop, c(1L, 0L, 0L, 1L))
})

test_that("notes record set and cleared transitions only", {
  out <- run_partial(dplyr::bind_rows(
    prow("cox1", "PCG", 200, 400, start_codon = "ATG", stop_codon = "TAA",
         partial_start = 1L, notes = "kept"),
    prow("cox1", "PCG", 401, 600, start_codon = "ATG", stop_codon = "AAA",
         partial_stop = NA_integer_),
    prow("cox1", "PCG", 601, 700, start_codon = "ATG", stop_codon = "TAA")
  ))
  expect_equal(out$notes[1], "kept; 5' partial cleared by curation")
  expect_equal(out$notes[2], "3' partial set by curation")
  expect_true(is.na(out$notes[3]))
  expect_equal(out$partial_start[1], 0L)
  expect_equal(out$partial_stop[2], 1L)
})

test_that("missing partial columns count as 0", {
  ann <- prow("trnF", "tRNA", 1, 70)
  ann$partial_start <- NULL
  ann$partial_stop <- NULL
  out <- run_partial(ann)
  expect_equal(out$partial_start, 1L)
  expect_equal(out$partial_stop, 0L)
  expect_equal(out$notes, "5' partial set by curation")
})

test_that("join_group_id parses JOIN markers", {
  expect_equal(
    join_group_id(c("JOIN: mode=exon group=3", "JOIN: mode=frameshift group=12 note=x", NA, "x; JOIN: mode=exon group=1")),
    c(3L, 12L, NA, NA)
  )
})

test_that("circular JOIN group: gap tie keeps the lower member 5'", {
  out <- run_partial(dplyr::bind_rows(
    prow("nad5", "PCG", 100, 300, contig = "C.1.1", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1"),
    prow("nad5", "PCG", 600, 800, contig = "C.1.1", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=exon group=1")
  ))
  expect_equal(out$partial_start, c(1L, 0L))
  expect_equal(out$partial_stop, c(0L, 1L))
})

test_that("circular JOIN group: overlapping members keep the lower member 5'", {
  out <- run_partial(dplyr::bind_rows(
    prow("nad5", "PCG", 100, 300, contig = "C.1.1", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=frameshift group=1"),
    prow("nad5", "PCG", 299, 500, contig = "C.1.1", start_codon = "CTG", stop_codon = "AAA",
         notes = "JOIN: mode=frameshift group=1")
  ))
  expect_equal(out$partial_start, c(1L, 0L))
  expect_equal(out$partial_stop, c(0L, 1L))
})

test_that("frame-offset PCG stays 5' partial even with a valid start codon", {
  ann <- dplyr::bind_rows(
    prow("cox1", "PCG", 200, 400, start_codon = "ATG", stop_codon = "TAA", partial_start = 1L),
    prow("nad5", "PCG", 401, 600, start_codon = "ATG", stop_codon = "TAA"),
    prow("nad6", "PCG", 601, 800, start_codon = "ATG", stop_codon = "TAA")
  )
  ann$codon_start <- c(2L, 3L, 1L)
  out <- run_partial(ann)
  expect_equal(out$partial_start, c(1L, 1L, 0L))
  expect_false(any(grepl("cleared", out$notes)))
  expect_equal(out$notes[2], "5' partial set by curation")
})

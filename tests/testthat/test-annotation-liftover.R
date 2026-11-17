set.seed(42)
rnd <- function(n) paste(sample(c("A", "C", "G", "T"), n, TRUE), collapse = "")
mutate <- function(s, sub = 0.08, indel = 0.02) {
  x <- strsplit(s, "")[[1]]
  mid <- 4:(length(x) - 3)
  i <- mid[runif(length(mid)) < sub]
  x[i] <- sample(c("A", "C", "G", "T"), length(i), TRUE)
  d <- mid[runif(length(mid)) < indel / 2]
  ins <- setdiff(mid[runif(length(mid)) < indel / 2], d)
  for (p in sort(ins, decreasing = TRUE)) x <- append(x, sample(c("A", "C", "G", "T"), 1), p)
  if (length(d)) x <- x[-d]
  paste(x, collapse = "")
}
rc <- function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
asm <- rnd(5000)
L <- nchar(asm)
all_hits <- list()
lift <- function(...) {
  h <- liftover_map(...)
  all_hits[[length(all_hits) + 1]] <<- h
  h
}
expect_hit <- function(h, p1, p2, dir = "+") {
  expect_equal(h$pos1[1], p1)
  expect_equal(h$pos2[1], p2)
  expect_equal(h$direction[1], dir)
}

test_that("liftover_map recovers forward, reverse, long, and wrapping genes", {
  h <- lift(mutate(substr(asm, 1001, 1070), indel = 0), asm, TRUE)
  expect_hit(h, 1001, 1070)

  h <- lift(rc(mutate(substr(asm, 2001, 2070), indel = 0)), asm, TRUE)
  expect_hit(h, 2001, 2070, "-")

  h <- lift(mutate(substr(asm, 3001, 3900)), asm, TRUE)
  expect_hit(h, 3001, 3900)

  h <- lift(mutate(substr(paste0(asm, asm), 4971, 5040), indel = 0), asm, TRUE)
  expect_hit(h, 4971, 40)
  expect_gt(h$pos1[1], h$pos2[1])
})

test_that("liftover_map returns both copies of a duplicated gene, best first", {
  dup <- rnd(70)
  a2 <- asm
  substr(a2, 501, 570) <- mutate(dup, 0.03, 0)
  substr(a2, 2501, 2570) <- mutate(dup, 0.12, 0)
  h <- lift(dup, a2, TRUE)
  expect_equal(nrow(h), 2)
  expect_equal(h$pos1, c(501, 2501))
  expect_equal(h$pos2, c(570, 2570))
  expect_gt(h$score[1], h$score[2])
})

test_that("liftover_map finds nothing for absent genes", {
  expect_equal(nrow(lift(rnd(70), asm, TRUE)), 0)
  expect_equal(nrow(lift(rnd(900), asm, TRUE)), 0)
})

test_that("liftover_map on linear assemblies", {
  h <- lift(mutate(substr(asm, 4931, 5000), indel = 0), asm, FALSE)
  expect_hit(h, 4931, 5000)
  h <- lift(rnd(6000), asm, FALSE)
  expect_equal(nrow(h), 0)
  expect_named(h, c("pos1", "pos2", "direction", "pctid", "gap_pct", "score", "aln_gene", "aln_sample"))
})

test_that("liftover_map tolerates N runs and ambiguity codes", {
  a3 <- asm
  substr(a3, 1501, 1550) <- strrep("N", 50)
  for (p in c(100, 900, 2533, 4000)) substr(a3, p, p) <- "R"
  for (p in c(300, 2545, 4500)) substr(a3, p, p) <- "Y"
  substr(a3, 4800, 4800) <- "?"
  h <- lift(mutate(substr(asm, 2501, 2570), indel = 0), a3, TRUE)
  expect_hit(h, 2501, 2570)
})

test_that("every returned hit passes the floors", {
  h <- do.call(rbind, all_hits)
  expect_gt(nrow(h), 0)
  expect_true(all(h$pctid >= 80))
  expect_true(all(h$gap_pct <= 10))
})

test_that("liftover_ref_gene_seq slices in gene orientation", {
  s <- "AAACCCGGGTTA"
  expect_equal(liftover_ref_gene_seq(s, 4, 6, "+"), "CCC")
  expect_equal(liftover_ref_gene_seq(s, 4, 9, "-"), "CCCGGG")
  expect_equal(liftover_ref_gene_seq(s, 7, 12, "-"), "TAACCC")
  expect_equal(liftover_ref_gene_seq(s, 11, 2, "+"), "TAAA")
})

test_that("liftover_blocked_reason", {
  br <- data.frame(
    gene = c("ctrl", "nad5", "nad5", "trnL", "trnL", "cox1"),
    type = c("ctrl", "PCG", "PCG", "tRNA", "tRNA", "PCG"),
    pos1 = c(1, 100, 900, 2000, 3000, 4000)
  )
  expect_type(liftover_blocked_reason(br, "ctrl", 1), "character")
  expect_false(is.na(liftover_blocked_reason(br, "ctrl", 1)))
  expect_false(is.na(liftover_blocked_reason(br, "nad5", 900)))
  expect_true(is.na(liftover_blocked_reason(br, "trnL", 3000)))
  expect_true(is.na(liftover_blocked_reason(br, "cox1", 4000)))
})

test_that("liftover_overlaps", {
  ann <- data.frame(
    gene = c("trnA", "trnB", "trnW", "cox1_DELETED_"),
    pos1 = c(171, 196, 4981, 1),
    pos2 = c(240, 265, 50, 400)
  )
  expect_equal(nrow(liftover_overlaps(1000, 1069, ann, 5000)), 0)
  o <- liftover_overlaps(106, 175, ann, 5000)
  expect_equal(o$gene, "trnA")
  expect_equal(o$bp, 5)
  expect_false(o$substantial)
  o <- liftover_overlaps(127, 196, ann, 5000)
  expect_equal(o$gene, c("trnA", "trnB"))
  expect_equal(o$bp, c(26, 1))
  expect_equal(o$substantial, c(TRUE, FALSE))
  o <- liftover_overlaps(21, 90, ann, 5000)
  expect_equal(o$gene, "trnW")
  expect_equal(o$bp, 30)
  expect_true(o$substantial)
})

test_that("liftover_clamp", {
  expect_equal(liftover_clamp(50, 40), list(min_pctid = 80, max_gap_pct = 10))
  expect_equal(liftover_clamp(95, 2), list(min_pctid = 95, max_gap_pct = 2))
})

test_that("liftover_map returns the scored alignment in gene orientation", {
  g1 <- mutate(substr(asm, 3001, 3400))
  g2 <- rc(mutate(substr(asm, 2001, 2070), indel = 0))
  g3 <- mutate(substr(paste0(asm, asm), 4971, 5040), indel = 0)
  for (g in list(g1, g2, g3)) {
    h <- liftover_map(g, asm, TRUE)
    expect_true(all(c("aln_gene", "aln_sample") %in% names(h)))
    expect_equal(nchar(h$aln_gene[1]), nchar(h$aln_sample[1]))
    expect_equal(gsub("-", "", h$aln_gene[1]), g)
  }
  h <- liftover_map(g2, asm, TRUE)
  expect_equal(h$direction[1], "-")
  expect_equal(gsub("-", "", h$aln_sample[1]), rc(substr(asm, 2001, 2070)))
})

test_that("identity counts only ACGT matches", {
  a4 <- asm
  substr(a4, 1031, 1040) <- strrep("N", 10)
  g <- substr(a4, 1001, 1070)
  h <- liftover_map(g, a4, TRUE)
  expect_hit(h, 1001, 1070)
  expect_lt(h$pctid[1], 100)
  expect_equal(h$pctid[1], 100 * 60 / 70)
})

test_that("a copy inside the appended head is found once, with the other copy", {
  dup <- rnd(70)
  a5 <- asm
  substr(a5, 1, 70) <- mutate(dup, 0.03, 0)
  substr(a5, 2501, 2570) <- mutate(dup, 0.12, 0)
  h <- liftover_map(dup, a5, TRUE)
  expect_equal(sort(h$pos1), c(1, 2501))
  expect_equal(sort(h$pos2), c(70, 2570))
})

test_that("liftover_overlaps: wrapping hit against non-wrapping annotations", {
  ann <- data.frame(gene = c("trnX", "trnY", "trnZ"), pos1 = c(10, 4950, 100), pos2 = c(60, 4995, 200))
  o <- liftover_overlaps(4990, 30, ann, 5000)
  expect_equal(o$gene, c("trnX", "trnY"))
  expect_equal(o$bp, c(21, 6))
  expect_equal(o$substantial, c(TRUE, FALSE))
})

test_that("liftover_ref_gene_seq wraps an end stored past the reference length", {
  s <- "AAACCCGGGTTA"
  expect_equal(liftover_ref_gene_seq(s, 11, 14, "+"), "TAAA")
  expect_equal(liftover_ref_gene_seq(s, 11, 14, "-"), "TTTA")
  ref <- strrep("ACGTTGCA", 2067)
  expect_equal(nchar(liftover_ref_gene_seq(ref, 16165, 16686, "-")), 522)
})

# PCG built from sense codons of the vertebrate mito code (2)
codons <- setdiff(
  apply(expand.grid(c("A", "C", "G", "T"), c("A", "C", "G", "T"), c("A", "C", "G", "T")), 1, paste, collapse = ""),
  c("TAA", "TAG", "AGA", "AGG")
)
pcg <- paste0("ATG", paste(sample(codons, 98, TRUE), collapse = ""), "TAA")
asm_p <- paste0(substr(asm, 1, 1000), pcg, substr(asm, 1301, 5000))
labs <- c("NC_1 nad1", "S1 1,001-1,300 (+)")
aln_for <- function(gene, a, starts = "ATG", circular = TRUE) {
  h <- liftover_map(gene, a, circular)
  ref <- data.frame(gene = "nad1", type = "PCG", pos1 = 1001, pos2 = 1300, direction = "+")
  liftover_alignment(h[1, ], ref, asm_p, a, "PCG", starts, 2L, 2L, labs)
}

test_that("liftover_alignment: tRNA gives the two-row nucleotide alignment", {
  h <- liftover_map(mutate(substr(asm, 1001, 1070), indel = 0), asm, TRUE)
  ref <- data.frame(gene = "trnF", type = "tRNA", pos1 = 1001, pos2 = 1070, direction = "+")
  x <- liftover_alignment(h[1, ], ref, asm, asm, "tRNA", "ATG", 2L, 2L, c("a", "b"))
  expect_equal(x$kind, "nt")
  expect_length(x$aln, 2)
  expect_equal(names(x$aln), c("a", "b"))
  expect_equal(as.character(x$aln), c(a = h$aln_gene[1], b = h$aln_sample[1]))
  expect_length(x$notes, 0)
})

test_that("liftover_alignment: PCG with a valid start gives a protein alignment", {
  x <- aln_for(pcg, asm_p)
  expect_equal(x$kind, "aa")
  expect_s4_class(x$aln, "AAStringSet")
  expect_equal(names(x$aln), labs)
  expect_equal(gsub("-", "", as.character(x$aln[[2]])), gsub("-", "", as.character(x$aln[[1]])))
  expect_length(x$notes, 0)
})

test_that("liftover_alignment: PCG without a valid start falls back to nucleotides", {
  a <- asm_p
  substr(a, 1001, 1003) <- "CCC"
  x <- aln_for(pcg, a)
  expect_equal(x$kind, "nt")
  expect_s4_class(x$aln, "DNAStringSet")
  expect_match(x$notes, "does not start with a valid start codon \\(CCC\\)", all = FALSE)
  x <- aln_for(pcg, asm_p, starts = "GTG")
  expect_equal(x$kind, "nt")
})

test_that("liftover_alignment: unknown start codons show nucleotides with a note", {
  for (sc in list(NULL, character())) {
    x <- aln_for(pcg, asm_p, starts = sc)
    expect_equal(x$kind, "nt")
    expect_match(x$notes, "Start codons for this gene are not known", all = FALSE)
    expect_false(any(grepl("does not start", x$notes)))
  }
})

test_that("liftover_alignment: internal stops are counted", {
  a <- asm_p
  substr(a, 1151, 1153) <- "TAA"
  x <- aln_for(pcg, a)
  expect_equal(x$kind, "aa")
  expect_match(x$notes, "^1 internal stop codon in the mapped region\\.$", all = FALSE)
})

test_that("liftover_alignment: wrapping and minus-strand placements", {
  aw <- paste0(substr(pcg, 151, 300), substr(asm, 1, 4700), substr(pcg, 1, 150))
  h <- liftover_map(pcg, aw, TRUE)
  expect_hit(h, 4851, 150)
  x <- aln_for(pcg, aw)
  expect_equal(x$kind, "aa")
  expect_length(x$notes, 0)
  am <- paste0(substr(asm, 1, 1000), rc(pcg), substr(asm, 1301, 5000))
  h <- liftover_map(pcg, am, TRUE)
  expect_equal(h$direction[1], "-")
  x <- aln_for(pcg, am)
  expect_equal(x$kind, "aa")
  expect_length(x$notes, 0)
})

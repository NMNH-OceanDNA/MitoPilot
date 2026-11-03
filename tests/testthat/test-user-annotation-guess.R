test_that("guess_gene maps names", {
  g <- function(n, t, p = NA) guess_gene(n, t, p)
  chk <- function(n, t, p, gene, type) {
    r <- guess_gene(n, t, p)
    expect_equal(r$gene, gene)
    expect_equal(r$type, type)
  }
  chk("COX1", "CDS", NA, "cox1", "PCG")
  chk("YP_000001.1", "CDS", "cytochrome c oxidase subunit I", "cox1", "PCG")
  chk("COI", "CDS", NA, "cox1", "PCG")
  chk("ATPase 6", "CDS", NA, "atp6", "PCG")
  chk("ND4L", "gene", NA, "nad4l", "PCG")
  chk("16S", "rRNA", NA, "rrnL", "rRNA")
  chk("RNR1", "rRNA", "12S ribosomal RNA", "rrnS", "rRNA")
  chk("tRNA-Leu (CUN)", "tRNA", NA, "trnL", "tRNA")
  chk("trnL2(taa)", "tRNA", NA, "trnL", "tRNA")
  chk("tRNA-Ser(uca)", NA, NA, "trnS", "tRNA")
  chk("TRNF", "gene", "tRNA-Phe", "trnF", "tRNA")
  chk("cox1", NA, NA, "cox1", "PCG")
  chk("D-loop", "misc_feature", NA, NA_character_, NA_character_)
  chk("OH", "gene", NA, NA_character_, NA_character_)
  chk("ORF1", "CDS", NA, NA_character_, NA_character_)
})

test_that("guess_gene is vectorised", {
  r <- guess_gene(c("COX1", "OH"), c("CDS", "gene"))
  expect_equal(nrow(r), 2L)
  expect_equal(r$gene, c("cox1", NA))
})

test_that("parse_anticodon", {
  expect_equal(parse_anticodon("tRNA-Ser (uca)"), "TCA")
  expect_equal(parse_anticodon("anticodon:gaa"), "GAA")
  expect_equal(parse_anticodon("(pos:1234..1236,aa:Phe,seq:gaa)"), "GAA")
  expect_equal(parse_anticodon("tRNA-Leu-UAA"), "TAA")
  expect_true(is.na(parse_anticodon("tRNA-Phe")))
  expect_true(is.na(parse_anticodon(NA)))
})

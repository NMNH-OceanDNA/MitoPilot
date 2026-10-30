# tests/testthat/test-synteny-payload.R
syn_ann <- function() data.frame(
  type = c("PCG", "tRNA", "PCG"), gene = c("nad1", "trnF", "cox1_DELETED_1"),
  pos1 = c(10L, 40L, 0L), pos2 = c(39L, 108L, 0L), direction = c("+", "-", "+"),
  partial_start = c(1L, 0L, 0L), partial_stop = c(0L, NA, 0L),
  stringsAsFactors = FALSE)
syn_ref <- function() data.frame(
  gene = c("nad1", "trnF"), type = c("PCG", "tRNA"), pos1 = c(5L, 90L), pos2 = c(34L, 20L),
  direction = c("+", "+"), ref_length = 100L, partial_start = c(0L, 1L), partial_stop = c(1L, 0L),
  stringsAsFactors = FALSE)
syn_aln <- function(rc = FALSE) data.frame(
  aligned_sample = "acg-t", aligned_ref = "ACGTT", rotation = 4L, ref_length = 100L,
  ref_rc = rc, stringsAsFactors = FALSE)

test_that("payload carries alignment, labels, and live features", {
  p <- synteny_payload(syn_ann(), syn_ref(), syn_aln(TRUE), "S1", "NC_1", "Fish sp.",
                       selected = 2L, version = 3L)
  expect_equal(p[c("sampleLabel", "refAcc", "refLabel")],
               list(sampleLabel = "S1", refAcc = "NC_1", refLabel = "Fish sp."))
  expect_true(p$refRc)
  expect_equal(p$s, "ACG-T"); expect_equal(p$r, "ACGTT")
  expect_equal(p$refLength, 100L); expect_equal(p$rotation, 4L)
  expect_equal(p$version, 3L); expect_equal(p$selected, 2L)
  expect_equal(vapply(p$sampleFeatures, `[[`, integer(1), "row"), c(1L, 2L))
  expect_equal(p$sampleFeatures[[2]][c("gene", "type", "pos1", "pos2", "dir")],
               list(gene = "trnF", type = "tRNA", pos1 = 40L, pos2 = 108L, dir = "-"))
  expect_equal(p$refFeatures[[2]][c("gene", "pos1", "pos2", "dir")],
               list(gene = "trnF", pos1 = 90L, pos2 = 20L, dir = "+"))
})

test_that("partial flags ride along as partial5/partial3, NA and missing columns read as complete", {
  p <- synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x")
  expect_equal(lapply(p$sampleFeatures, `[`, c("partial5", "partial3")),
               list(list(partial5 = TRUE, partial3 = FALSE), list(partial5 = FALSE, partial3 = FALSE)))
  expect_equal(lapply(p$refFeatures, `[`, c("partial5", "partial3")),
               list(list(partial5 = FALSE, partial3 = TRUE), list(partial5 = TRUE, partial3 = FALSE)))
  old_ref <- syn_ref()[, c("gene", "type", "pos1", "pos2", "direction", "ref_length")]
  p2 <- synteny_payload(syn_ann(), old_ref, syn_aln(), "S1", "NC_1", "x")
  expect_equal(p2$refFeatures[[1]][c("partial5", "partial3")], list(partial5 = FALSE, partial3 = FALSE))
})

test_that("selection outside the live rows is dropped and missing alignment gives NULL", {
  expect_null(synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x", selected = 3L)$selected)
  expect_false(synteny_payload(syn_ann(), syn_ref(), syn_aln(), "S1", "NC_1", "x")$refRc)
  expect_null(synteny_payload(syn_ann(), syn_ref(), NULL, "S1", "NC_1", "x"))
  expect_null(synteny_payload(syn_ann(), syn_ref(), syn_aln()[0, ], "S1", "NC_1", "x"))
  empty <- syn_aln(); empty$aligned_sample <- ""
  expect_null(synteny_payload(syn_ann(), syn_ref(), empty, "S1", "NC_1", "x"))
  expect_length(synteny_payload(syn_ann(), syn_ref()[0, ], syn_aln(), "S1", "NC_1", "x")$refFeatures, 0)
})

test_that("viewer UI carries the namespaced canvas, tip, and the four controls", {
  html <- as.character(synteny_viewer_ui(shiny::NS("a")))
  expect_match(html, 'id="a-synteny_canvas"', fixed = TRUE)
  expect_match(html, 'id="a-synteny_tip"', fixed = TRUE)
  for (a in c("fit", "whole", "zoom_in", "zoom_out")) {
    expect_match(html, sprintf('data-mpsyn="%s"', a), fixed = TRUE)
  }
})

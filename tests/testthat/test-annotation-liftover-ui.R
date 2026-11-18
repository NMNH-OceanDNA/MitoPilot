ann <- data.frame(
  gene = c("trnV", "rrnS", "trnF_DELETED_1"),
  pos1 = c(1190L, 3000L, 1200L),
  pos2 = c(1260L, 3950L, 1300L)
)
hits <- data.frame(
  pos1 = c(1204L, 9877L, 4990L), pos2 = c(2153L, 10820L, 30L),
  direction = c("+", "+", "-"),
  pctid = c(93.14, 84.96, 81.5), gap_pct = c(1.21, 4.1, 7.9), score = c(900, 700, 600)
)

test_that("hit table formats cells, flags overlaps, and applies thresholds", {
  tb <- liftover_hit_table(hits, ann, 16000, 85, 5)
  expect_equal(tb$position, c("1,204-2,153", "9,877-10,820", "4,990-30"))
  expect_equal(tb$strand, c("+", "+", "-"))
  expect_equal(tb$identity, c("93.1%", "84.9%", "81.5%"))
  expect_equal(tb$gaps, c("1.3%", "4.1%", "7.9%"))
  expect_equal(tb$overlaps, c("trnV (57 bp)", "none", "none"))
  expect_equal(tb$clash, c("trnV", "", ""))
  expect_equal(tb$pass, c(TRUE, FALSE, FALSE))
})

test_that("hit table clamps thresholds to the hard floors", {
  expect_equal(liftover_hit_table(hits, ann, 16000, 50, 50)$pass, c(TRUE, TRUE, TRUE))
  expect_equal(liftover_hit_table(hits, ann, 16000, 84, 5)$pass, c(TRUE, TRUE, FALSE))
})

test_that("hit table handles zero hits", {
  tb <- liftover_hit_table(hits[0, ], ann, 16000, 85, 5)
  expect_equal(nrow(tb), 0L)
})

test_that("confirm text names clashing genes with Oxford commas", {
  expect_null(liftover_confirm_text("trnL", character(), FALSE))
  expect_equal(
    liftover_confirm_text("trnL", "trnV", FALSE),
    "This placement overlaps trnV by more than a quarter of either gene's length. Add trnL anyway?"
  )
  expect_match(liftover_confirm_text("trnL", c("a", "b", "c"), TRUE),
               "overlaps a, b, and c by .* already has an annotation named trnL\\. Add trnL anyway\\?$")
  expect_match(liftover_confirm_text("trnL", c("a", "b"), FALSE), "overlaps a and b by", fixed = FALSE)
})

test_that("is_liftover_row only matches liftover tool tags", {
  expect_true(is_liftover_row("NC_012345.1_liftover"))
  expect_false(is_liftover_row("MITOS2"))
  expect_false(is_liftover_row(NA_character_))
  expect_false(is_liftover_row(character(0)))
})

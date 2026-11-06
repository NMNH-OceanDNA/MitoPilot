# tests/testthat/test-user-annot-track.R
rot <- function(s, n) paste0(substr(s, n + 1, nchar(s)), substr(s, 1, n))
set.seed(1)
S <- paste(sample(c("A", "C", "G", "T"), 2000, TRUE), collapse = "")

test_that("find_frame: identity, rotation, rc, probe over the new origin", {
  expect_equal(find_frame(S, S), list(offset = 0L, rc = FALSE))
  expect_equal(find_frame(S, rot(S, 500)), list(offset = -500L, rc = FALSE))
  expect_equal(find_frame(S, revcomp_chr(S)), list(offset = 0L, rc = TRUE))
  expect_equal(find_frame(S, rot(S, 1990)), list(offset = -1990L, rc = FALSE))
  expect_equal(find_frame(S, rot(S, 30)), list(offset = -30L, rc = FALSE))
  expect_equal(find_frame(S, revcomp_chr(rot(S, 500))), list(offset = -500L, rc = TRUE))
})

test_that("find_frame: trimmed start uses a later probe; repeats and unrelated give NULL", {
  expect_equal(find_frame(S, substr(S, 101, 1900)), list(offset = -100L, rc = FALSE))
  expect_null(find_frame(S, paste0(S, S)))
  expect_null(find_frame(S, paste(rev(strsplit(S, "")[[1]]), collapse = "")))
})

test_that("map_feature: rotation keeps strand, rc flips, wrap allowed", {
  f <- list(offset = -500L, rc = FALSE)
  expect_equal(map_feature(600L, 700L, "+", f, 2000L, 2000L)[c("pos1", "pos2", "dir")],
               list(pos1 = 100L, pos2 = 200L, dir = "+"))
  expect_equal(map_feature(400L, 600L, "+", f, 2000L, 2000L)[c("pos1", "pos2")],
               list(pos1 = 1900L, pos2 = 100L))
  r <- list(offset = 0L, rc = TRUE)
  expect_equal(map_feature(10L, 30L, "+", r, 2000L, 2000L)[c("pos1", "pos2", "dir")],
               list(pos1 = 1971L, pos2 = 1991L, dir = "-"))
  expect_false(map_feature(10L, 30L, "+", r, 2000L, 2000L)$clipped)
})

test_that("map_feature: linear shift clips, flags, and drops", {
  f <- list(offset = -100L, rc = FALSE)
  expect_equal(map_feature(50L, 150L, "+", f, 2000L, 1800L),
               list(pos1 = 1L, pos2 = 50L, dir = "+", clipped = TRUE))
  expect_false(map_feature(500L, 600L, "+", f, 2000L, 1800L)$clipped)
  expect_null(map_feature(10L, 90L, "+", f, 2000L, 1800L))
  expect_equal(map_feature(1950L, 300L, "+", f, 2000L, 1800L)[c("pos1", "pos2", "clipped")],
               list(pos1 = 1L, pos2 = 200L, clipped = TRUE))
})

test_that("truncated_match separates truncated calls from dissimilar ones", {
  h <- function(pctid, gl = 0L, gt = 0L) data.frame(pctid = pctid, gap_leading = gl, gap_trailing = gt)
  expect_true(truncated_match(h(100, gl = 76L), 90))
  expect_true(truncated_match(h(c(50, 95), gt = c(0L, 12L)), 90))
  expect_false(truncated_match(h(60, gl = 76L), 90))
  expect_false(truncated_match(h(100), 90))
  expect_false(truncated_match(h(100, gl = -4L), 90))
  expect_false(truncated_match(h(numeric(0), integer(0), integer(0)), 90))
})

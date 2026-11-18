# update_counts() tallies: a JOIN group is one gene copy.

join_rows <- function() {
  data.frame(
    contig = c("c1", "c1", "c1", "c2", "c2"),
    type = c("PCG", "PCG", "tRNA", "PCG", "PCG"),
    gene = c("nad5", "nad5", "trnF", "nad5", "nad5"),
    notes = c("JOIN: mode=exon group=1", "JOIN: mode=exon group=1", NA,
              "JOIN: mode=exon group=1", NA),
    stringsAsFactors = FALSE
  )
}

test_that("count_feature_types counts a JOIN group once", {
  a <- join_rows()[1:3, ]
  cnt <- count_feature_types(a)
  expect_equal(cnt[["PCG"]], 1L)
  expect_equal(cnt[["tRNA"]], 1L)
  expect_equal(cnt[["ORF"]], 0L)
  # same group id on another contig, and a non-member, count separately
  expect_equal(count_feature_types(join_rows())[["PCG"]], 3L)
})

test_that("compute_missing_extra counts a JOIN group once", {
  rules <- list(nad5 = list(type = "PCG", count = 1))
  me <- compute_missing_extra(join_rows()[1:3, ], rules, list(PCG = list()))
  expect_true(is.na(me$extra))
  expect_true(is.na(me$missing))
  me <- compute_missing_extra(join_rows(), rules, list(PCG = list()))
  expect_equal(me$extra, "nad5")
})

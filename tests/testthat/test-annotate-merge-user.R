test_that("same-gene overlap drops the tool row", {
  user <- data.frame(contig = "c", gene = "nad1", pos1 = 100L, pos2 = 1000L)
  tool <- data.frame(contig = "c", gene = "nad1", pos1 = 60L, pos2 = 1000L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 0L)
})
test_that("small edge overlap of a different gene is kept", {
  user <- data.frame(contig = "c", gene = "nad1", pos1 = 100L, pos2 = 1000L)
  tool <- data.frame(contig = "c", gene = "trnA", pos1 = 996L, pos2 = 1065L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 1L)
})
test_that("over 10 percent overlap of a different gene is dropped", {
  user <- data.frame(contig = "c", gene = "nad1", pos1 = 100L, pos2 = 1000L)
  tool <- data.frame(contig = "c", gene = "nad2", pos1 = 900L, pos2 = 1500L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 0L)
})
test_that("origin-spanning tool row is measured wrapped", {
  user <- data.frame(contig = "c", gene = "cob", pos1 = 1L, pos2 = 300L)
  tool <- data.frame(contig = "c", gene = "nad6", pos1 = 4900L, pos2 = 150L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 0L)
})
test_that("MITOS OH counts as user ctrl", {
  user <- data.frame(contig = "c", gene = "ctrl", pos1 = 2000L, pos2 = 2900L)
  tool <- data.frame(contig = "c", gene = "OH", pos1 = 2850L, pos2 = 3800L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 0L)
})
test_that("same gene elsewhere on the contig is kept, other contigs ignored", {
  user <- data.frame(contig = "c", gene = "cox1", pos1 = 100L, pos2 = 1600L)
  tool <- data.frame(contig = c("c", "d"), gene = c("cox1", "cox1"),
                     pos1 = c(3000L, 100L), pos2 = c(3400L, 1600L))
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L, d = 5000L))), 2L)
})
test_that("empty inputs are safe", {
  user <- data.frame(contig = character(), gene = character(), pos1 = integer(), pos2 = integer())
  tool <- data.frame(contig = "c", gene = "trnA", pos1 = 1L, pos2 = 70L)
  expect_equal(nrow(drop_overlapping_tool_rows(tool, user, c(c = 5000L))), 1L)
  expect_equal(nrow(drop_overlapping_tool_rows(user, tool, c(c = 5000L))), 0L)
})

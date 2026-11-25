test_that("install_cli writes an executable wrapper bound to an Rscript", {
  d <- withr::local_tempdir()
  p <- suppressMessages(install_cli(dir = d, rscript = "/opt/R/bin/Rscript"))
  expect_true(file.exists(p))
  expect_equal(basename(p), "mitopilot")
  expect_true(any(grepl('exec "/opt/R/bin/Rscript"', readLines(p), fixed = TRUE)))
  expect_equal(unname(file.access(p, 1)), 0)
})

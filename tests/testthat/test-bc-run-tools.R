bc_run_tools_db <- function(td, user_annot = TRUE) {
  suppressMessages(new_test_project_userAnnot(
    path = td, no_raw_data = TRUE, Rproj = FALSE, fetch_metadata = FALSE))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  DBI::dbExecute(con, "ALTER TABLE annotate_opts DROP COLUMN run_tools")
  DBI::dbExecute(con, "UPDATE annotate_opts SET use_mitos = 0, use_trnaScan = 0")
  if (!user_annot) DBI::dbExecute(con, "DROP TABLE gene_name_map")
  DBI::dbDisconnect(con)
  suppressMessages(backwards_compatibility(path = td, update_config = FALSE))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbReadTable(con, "annotate_opts")
}

test_that("migration adds run_tools off and restores tool defaults on user-annotation projects", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  ao <- bc_run_tools_db(withr::local_tempdir())
  expect_equal(ao$run_tools, 0L)
  expect_equal(ao$use_mitos, 1L)
  expect_equal(ao$use_trnaScan, 1L)
})

test_that("migration adds run_tools on for other projects", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  ao <- bc_run_tools_db(withr::local_tempdir(), user_annot = FALSE)
  expect_equal(ao$run_tools, 1L)
  expect_equal(ao$use_mitos, 0L)
  expect_equal(ao$use_trnaScan, 0L)
})

test_that("new projects carry settings, review_log, and export option defaults", {
  con <- local_test_db()
  expect_equal(project_settings(con), list(min_mean_depth = 10, accept_linear = FALSE))
  expect_true(DBI::dbExistsTable(con, "review_log"))
  eo <- DBI::dbGetQuery(con, "SELECT * FROM export_opts WHERE export_opts = 'default'")
  expect_equal(c(eo$start_aa, eo$stop_aa, eo$ident_pct, eo$aa_alignments, eo$gene_export), c(10, 10, 60, 1, 0))
  expect_false(any(grepl("command-line", schema_gaps(con))))
})

test_that("migration adds the new pieces to an old database", {
  con <- local_test_db()
  DBI::dbExecute(con, "DROP TABLE project_settings")
  DBI::dbExecute(con, "DROP TABLE review_log")
  DBI::dbExecute(con, "ALTER TABLE export_opts DROP COLUMN ident_pct")
  expect_true(any(grepl("command-line", schema_gaps(con))))
  DBI::dbExecute(con, "INSERT INTO export_opts (export_opts, fasta_header) VALUES ('mine', 'x')")
  migrate_cli_schema(con)
  expect_equal(project_settings(con)$min_mean_depth, 10)
  expect_true("ident_pct" %in% DBI::dbListFields(con, "export_opts"))
  expect_equal(DBI::dbGetQuery(con, "SELECT ident_pct FROM export_opts WHERE export_opts = 'mine'")$ident_pct, 60)
  expect_false(any(grepl("command-line", schema_gaps(con))))
})

test_that("open_project fails cleanly outside a project", {
  expect_error(open_project(withr::local_tempdir()), class = "mitopilot_cli_error")
})

test_that("open_project returns a handle and migrates an old database", {
  td <- withr::local_tempdir()
  con <- local_test_db()
  db <- DBI::dbGetInfo(con)$dbname
  DBI::dbExecute(con, "DROP TABLE review_log")
  file.copy(db, file.path(td, ".sqlite"))
  writeLines(c("params {", "    asmbDir = 'NA'", "    publishDir = 'out'", "}"), file.path(td, ".config"))
  withr::local_options(MitoPilot.db = NULL)
  p <- suppressMessages(open_project(td))
  withr::defer(DBI::dbDisconnect(p$con))
  expect_false(p$user_asmb)
  expect_equal(p$dir_out, file.path(normalizePath(td), "out"))
  expect_true(DBI::dbExistsTable(p$con, "review_log"))
})

test_that("set_project_setting rejects unknown keys and stores known ones", {
  con <- local_test_db()
  expect_error(set_project_setting(con, "bogus", 1), class = "mitopilot_cli_error")
  set_project_setting(con, "accept_linear", TRUE)
  set_project_setting(con, "min_mean_depth", 25)
  expect_equal(project_settings(con), list(min_mean_depth = 25, accept_linear = TRUE))
})

test_that("run_active is NA with no runs and finds an unreported live run", {
  td <- withr::local_tempdir()
  expect_true(is.na(run_active(td)))
  dir.create(file.path(td, ".runs", "nextflow"), recursive = TRUE)
  dir.create(file.path(td, ".runs", "reports"), recursive = TRUE)
  log <- file.path(td, ".runs", "nextflow", "assemble_2026-01-01_10-00-00.nextflow.log")
  writeLines("Starting process > a", log)
  expect_equal(run_active(td), "assemble_2026-01-01_10-00-00")
  file.create(file.path(td, ".runs", "reports", "assemble_2026-01-01_10-00-00_report.dcf"))
  expect_true(is.na(run_active(td)))
  file.remove(file.path(td, ".runs", "reports", "assemble_2026-01-01_10-00-00_report.dcf"))
  Sys.setFileTime(log, Sys.time() - 7200)
  expect_true(is.na(run_active(td)))
})

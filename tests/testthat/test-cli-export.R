flagged_two <- function(...) list(
  flagged = data.frame(ID = c("B", "B"), gene = c("ND5", "COX1"), why = c("low_identity", "internal_stop")),
  result = NULL, report = "/p/out/export/g1/AA_alignments_g1.html")
flagged_none <- function(...) list(
  flagged = data.frame(ID = character(), gene = character(), why = character()),
  result = NULL, report = "/p/out/export/g1/AA_alignments_g1.html")

test_that("scripted export stops by default when samples are flagged", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) stop("must not export"))
  msgs <- testthat::capture_messages(st <- cli_cmd_export(c("g1", "--project", f$path)))
  expect_equal(st, 2L); expect_true(any(grepl("AA_alignments_g1.html", msgs, fixed = TRUE)))
})
test_that("--flagged drop exports only unflagged samples and exits 2", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B")); got <- NULL
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) got <<- list(...))
  expect_equal(suppressMessages(cli_cmd_export(c("g1", "--flagged", "drop", "--project", f$path))), 2L)
  expect_equal(got$IDs, "A"); expect_false(got$review)
})
test_that("--flagged keep exports all and logs an override per flagged unit", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) NULL)
  expect_equal(suppressMessages(cli_cmd_export(c("g1", "--flagged", "keep", "--yes", "--project", f$path))), 0L)
  lg <- DBI::dbGetQuery(f$con, "SELECT ID, stage, action, rules FROM review_log")
  expect_equal(c(lg$ID, lg$stage, lg$action), c("B", "export", "override"))
  expect_setequal(strsplit(lg$rules, ",")[[1]], c("low_identity", "internal_stop"))
})
test_that("ask without a terminal fails", {
  f <- local_export_fixture(group = "g1", ids = "A")
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("export", "g1", "--flagged", "ask", "--project", f$path))), 1L)
})

test_that("keep in a script without --yes exports and logs nothing", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) stop("must not export"))
  expect_equal(suppressMessages(cli_main(c("export", "g1", "--flagged", "keep", "--project", f$path))), 1L)
  expect_equal(DBI::dbGetQuery(f$con, "SELECT COUNT(*) n FROM review_log")$n, 0L)
})

test_that("ask menu: drop, keep, cancel", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  got <- NULL
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) got <<- list(...))
  ask <- function(input) {
    got <<- NULL
    withr::local_options(MitoPilot.cli_interactive = TRUE, MitoPilot.cli_input = textConnection(input))
    suppressMessages(cli_cmd_export(c("g1", "--project", f$path)))
  }
  expect_equal(ask("1"), 2L); expect_equal(got$IDs, "A")
  expect_equal(ask("3"), 2L); expect_null(got)
  expect_equal(ask(c("2", "y")), 0L); expect_null(got$IDs)
})

test_that("nothing flagged exports everything with the saved options", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B")); got <- NULL
  DBI::dbExecute(f$con, "UPDATE export_opts SET start_aa = 7, stop_aa = 8, ident_pct = 55, aa_alignments = 0, gene_export = 1")
  local_mocked_bindings(export_flags = flagged_none, export_files = function(...) got <<- list(...))
  msgs <- testthat::capture_messages(st <- cli_cmd_export(c("g1", "--project", f$path)))
  expect_equal(st, 0L)
  expect_null(got$IDs); expect_equal(got$group, "g1"); expect_equal(got$out_dir, f$dir_out)
  expect_equal(got[c("start_aa", "stop_aa", "ident_pct", "generateAAalignments", "gene_export", "review")],
               list(start_aa = 7, stop_aa = 8, ident_pct = 55, generateAAalignments = FALSE,
                    gene_export = TRUE, review = FALSE))
  expect_true(any(grepl("AA_alignments_g1.html", msgs, fixed = TRUE)))
})

test_that("every sample flagged with drop writes nothing and exits 2", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  two <- function(...) list(flagged = data.frame(ID = c("A", "B"), gene = "ND5", why = "low_identity"),
                            result = NULL, report = "r.html")
  local_mocked_bindings(export_flags = two, export_files = function(...) stop("must not export"))
  msgs <- testthat::capture_messages(st <- cli_cmd_export(c("g1", "--flagged", "drop", "--project", f$path)))
  expect_equal(st, 2L); expect_true(any(grepl("Nothing to export", msgs)))
})

test_that("blocks: unknown group, active run, invalid header, multiple paths", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  local_mocked_bindings(export_flags = flagged_none, export_files = function(...) stop("must not export"))
  run <- function(...) suppressMessages(cli_main(c("export", ..., "--project", f$path)))
  expect_equal(run("nope"), 1L)
  expect_equal(run("g1", "--flagged", "bogus"), 1L)
  local({
    local_mocked_bindings(run_active = function(path) "run_x")
    expect_equal(run("g1"), 1L)
  })
  DBI::dbExecute(f$con, "UPDATE export_opts SET fasta_header = '{seqid} {not_a_field}'")
  msg <- testthat::capture_messages(st <- cli_main(c("export", "g1", "--project", f$path)))
  expect_equal(st, 1L); expect_true(any(grepl("header template is invalid", msg)))

  m <- local_export_fixture(group = "g1", ids = "A",
                            annotate_units = data.frame(ID = "A", path = c(1, 2), scaffold = 1))
  msg <- testthat::capture_messages(st <- cli_main(c("export", "g1", "--project", m$path)))
  expect_equal(st, 1L); expect_true(any(grepl("more than one assembly path", msg)))
})

test_that("a fragmented sample needs --yes in a script", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"),
                            annotate_units = data.frame(ID = c("A", "A", "B"), path = 1, scaffold = c(1, 2, 1)))
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  got <- NULL
  local_mocked_bindings(export_flags = flagged_none, export_files = function(...) got <<- list(...))
  expect_equal(suppressMessages(cli_main(c("export", "g1", "--project", f$path))), 1L)
  expect_null(got)
  msgs <- testthat::capture_messages(st <- cli_main(c("export", "g1", "--yes", "--project", f$path)))
  expect_equal(st, 0L); expect_true(any(grepl("more than one assembly", msgs)))
  expect_equal(got$group, "g1")
})

test_that("unchecked sample-gene pairs are counted", {
  a <- data.frame(ID = c("A", "B", "C", "A", "B"), gene = c("nad5", "nad5", "nad5", "atp8", "cox1"),
                  translation = c("MK", NA, "MR", "MA", "MQ"))
  expect_equal(export_unchecked(a), 3L)
  expect_equal(export_unchecked(a[0, ]), 0L)
})

test_that("export options round-trip, with today's values as defaults", {
  con <- local_test_db()
  o <- get_export_opts(con)
  expect_equal(o[c("start_aa", "stop_aa", "ident_pct", "aa_alignments", "gene_export")],
               list(start_aa = 10, stop_aa = 10, ident_pct = 60, aa_alignments = TRUE, gene_export = FALSE))
  set_export_opts(con, start_aa = 5, stop_aa = 6, ident_pct = 70, aa_alignments = FALSE, gene_export = TRUE)
  o2 <- get_export_opts(con)
  expect_equal(o2[c("start_aa", "stop_aa", "ident_pct", "aa_alignments", "gene_export")],
               list(start_aa = 5, stop_aa = 6, ident_pct = 70, aa_alignments = FALSE, gene_export = TRUE))
  expect_equal(o2$fasta_header, o$fasta_header)
  set_export_opts(con, "{seqid} x {completeness}", "{seqid} g")
  expect_equal(get_export_opts(con)$start_aa, 5)
  DBI::dbExecute(con, "ALTER TABLE export_opts DROP COLUMN start_aa")
  expect_equal(get_export_opts(con)$start_aa, 10)
})

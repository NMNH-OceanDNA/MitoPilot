auto_q <- function(f, sql) DBI::dbGetQuery(f$con, sql)[[1]]

auto_fixture <- function(executor = "local", env = parent.frame()) {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE, depth = c(80, 2),
    annotate_units = data.frame(ID = c("A", "B"), path = 1, scaffold = 1),
    executor = executor, env = env)
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 1")
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 1")
  f
}

mock_runs <- function(f, calls) {
  function(project, workflow, ...) {
    calls$v <- c(calls$v, workflow)
    if (workflow == "assemble") {
      DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2")
    } else {
      DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 2, warnings = 0, problematic = 'no', partial = 'no'")
    }
    list(status = 0L, queued = 1L)
  }
}

test_that("auto accepts clean units, parks the rest, exports, and exits 2", {
  f <- auto_fixture()
  calls <- new.env(); exported <- NULL
  local_mocked_bindings(run_pipeline = mock_runs(f, calls),
    cli_cmd_export = function(args) { exported <<- args[1]; 0L })
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path, "--yes"))), 2L)
  expect_equal(calls$v, c("assemble", "annotate")); expect_equal(exported, "g1")
  st <- project_status(f$con)
  expect_equal(st$export_group[st$ID == "A"], "g1"); expect_true(is.na(st$export_group[st$ID == "B"]))
})

test_that("auto stops with exit 1 when the assemble run fails", {
  f <- auto_fixture()
  local_mocked_bindings(run_pipeline = function(...) list(status = 1L, queued = 1L))
  expect_equal(suppressMessages(cli_cmd_auto(c("--project", f$path, "--yes"))), 1L)
  expect_equal(auto_q(f, "SELECT COUNT(*) FROM export WHERE export_group IS NOT NULL"), 0)
})

test_that("a second auto launches no run and exports nothing new", {
  f <- auto_fixture()
  calls <- new.env(); n_exp <- 0L
  local_mocked_bindings(run_pipeline = mock_runs(f, calls),
    cli_cmd_export = function(args) {
      n_exp <<- n_exp + 1L
      DBI::dbExecute(f$con, "UPDATE export SET export_time_stamp = 1 WHERE export_group = ?", params = list(args[1]))
      0L
    })
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path))), 2L)
  calls$v <- character()
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path))), 2L)
  expect_length(calls$v, 0L); expect_equal(n_exp, 1L)
})

test_that("auto --dry-run launches and exports nothing", {
  f <- auto_fixture()
  local_mocked_bindings(run_pipeline = function(...) stop("must not launch"),
    cli_cmd_export = function(args) stop("must not export"))
  out <- capture.output(st <- suppressMessages(cli_cmd_auto(c("--dry-run", "--project", f$path))))
  expect_equal(st, 0L); expect_true(any(grepl("nextflow", out)))
  expect_equal(auto_q(f, "SELECT SUM(assemble_lock) FROM assemble"), 0)
})

test_that("auto --submit --dry-run prints a job script with auto ... --yes", {
  f <- auto_fixture("slurm")
  local_mocked_bindings(cli_self_command = function() "/x/mitopilot", run_pipeline = function(...) stop("must not launch"))
  out <- capture.output(st <- suppressMessages(cli_cmd_auto(c("--submit", "--dry-run", "--project", f$path))))
  ln <- grep("mitopilot auto", out, value = TRUE)
  expect_equal(st, 0L); expect_true(any(grepl("--yes", ln, fixed = TRUE)))
  expect_false(any(grepl("--submit", ln, fixed = TRUE)))
})

test_that("auto skips annotate when gene names are unconfirmed but still groups and exports", {
  f <- auto_fixture()
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 2, warnings = 0, problematic = 'no', partial = 'no'")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2")
  calls <- new.env(); exported <- character()
  local_mocked_bindings(run_pipeline = mock_runs(f, calls), annot_import_pending = function(...) TRUE,
    cli_cmd_export = function(args) { exported <<- c(exported, args[1]); 0L })
  suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path)))
  expect_length(calls$v, 0L); expect_equal(exported, "g1")
})

test_that("a failing export gives exit 2 overall", {
  f <- auto_fixture()
  calls <- new.env()
  local_mocked_bindings(run_pipeline = mock_runs(f, calls), cli_cmd_export = function(args) 2L)
  DBI::dbExecute(f$con, "UPDATE assemblies SET depth = '80 80 80 80 80'")
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path))), 2L)
  expect_equal(auto_q(f, "SELECT COUNT(*) FROM export WHERE export_group = 'g1'"), 2)
})

test_that("auto --flagged keep needs --yes in a script and stops before any run", {
  f <- auto_fixture()
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  local_mocked_bindings(run_pipeline = function(...) stop("must not launch"))
  expect_equal(suppressMessages(cli_main(c("auto", "--flagged", "keep", "--project", f$path))), 1L)
})

test_that("auto rejects a bad group name and a bad --flagged value", {
  f <- auto_fixture()
  expect_error(cli_cmd_auto(c("--group", "a b", "--project", f$path)), class = "mitopilot_cli_error")
  expect_error(cli_cmd_auto(c("--flagged", "ask", "--project", f$path)), class = "mitopilot_cli_error")
})

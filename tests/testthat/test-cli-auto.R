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
  expect_equal(st, 0L); expect_false(any(grepl("--yes", ln, fixed = TRUE)))
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

clean_fixture <- function(env = parent.frame()) {
  f <- auto_fixture(env = env)
  DBI::dbExecute(f$con, "UPDATE assemblies SET depth = '80 80 80 80 80'")
  f
}

test_that("an export that returns 2 gives exit 2 when nothing is parked, 0 otherwise", {
  f <- clean_fixture()
  calls <- new.env()
  local_mocked_bindings(run_pipeline = mock_runs(f, calls), cli_cmd_export = function(args) 2L)
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path))), 2L)
  expect_equal(auto_q(f, "SELECT COUNT(*) FROM export WHERE export_group = 'g1'"), 2)
  g <- clean_fixture()
  local_mocked_bindings(run_pipeline = mock_runs(g, calls), cli_cmd_export = function(args) 0L)
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", g$path))), 0L)
})

test_that("a failing annotate run exits 1 and groups nothing", {
  f <- clean_fixture()
  local_mocked_bindings(run_pipeline = function(project, workflow, ...) {
    if (workflow == "assemble") DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2")
    list(status = if (workflow == "annotate") 1L else 0L, queued = 1L)
  }, cli_cmd_export = function(args) stop("must not export"))
  expect_equal(suppressMessages(cli_cmd_auto(c("--project", f$path))), 1L)
  expect_equal(auto_q(f, "SELECT COUNT(*) FROM export WHERE export_group IS NOT NULL"), 0)
})

test_that("an export that errors gives exit 1 and later groups still export", {
  f <- clean_fixture()
  calls <- new.env(); done <- character()
  run <- mock_runs(f, calls)
  local_mocked_bindings(
    run_pipeline = function(project, workflow, ...) {
      r <- run(project, workflow, ...)
      if (workflow == "annotate") DBI::dbExecute(f$con, "UPDATE annotate SET partial = 'yes' WHERE ID = 'B'")
      r
    },
    cli_cmd_export = function(args) {
      if (args[1] == "g1-complete") cli_fail("bad header")
      done <<- c(done, args[1]); 0L
    })
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path))), 1L)
  expect_equal(done, "g1-partial")
})

test_that("an export that needs confirmation is held: exit 2, with the way forward", {
  f <- clean_fixture()
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  calls <- new.env()
  local_mocked_bindings(run_pipeline = mock_runs(f, calls),
    cli_cmd_export = function(args) { cli_confirm("Export anyway?", FALSE, detail = "NMNH values missing"); 0L })
  msgs <- testthat::capture_messages(st <- cli_cmd_auto(c("--group", "g1", "--project", f$path)))
  expect_equal(st, 2L)
  expect_true(any(grepl("NMNH values missing", msgs, fixed = TRUE)))
  expect_true(any(grepl("mitopilot export g1 --yes", msgs, fixed = TRUE)))
})

test_that("auto --yes is passed to export only when given", {
  f <- clean_fixture()
  calls <- new.env(); got <- list()
  local_mocked_bindings(run_pipeline = mock_runs(f, calls), cli_cmd_export = function(args) { got[[length(got) + 1L]] <<- args; 0L })
  suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path)))
  expect_false("--yes" %in% got[[1]])
  DBI::dbExecute(f$con, "UPDATE export SET export_group = NULL")
  suppressMessages(cli_cmd_auto(c("--group", "g1", "--yes", "--project", f$path)))
  expect_true("--yes" %in% got[[2]])
})

test_that("auto --flagged keep confirms once and tells export it is confirmed", {
  f <- clean_fixture()
  calls <- new.env(); got <- NULL
  withr::local_options(MitoPilot.cli_interactive = TRUE, MitoPilot.cli_input = textConnection("y"))
  local_mocked_bindings(run_pipeline = mock_runs(f, calls), cli_cmd_export = function(args) { got <<- args; 0L })
  suppressMessages(cli_cmd_auto(c("--flagged", "keep", "--group", "g1", "--project", f$path)))
  expect_true("--flagged-confirmed" %in% got); expect_false("--yes" %in% got)
})

test_that("auto --submit adds --yes to the job only when given", {
  f <- auto_fixture("slurm")
  local_mocked_bindings(cli_self_command = function() "/x/mitopilot")
  line <- function(...) {
    out <- capture.output(suppressMessages(cli_cmd_auto(c("--submit", "--dry-run", ..., "--project", f$path))))
    grep("mitopilot auto", out, value = TRUE)
  }
  expect_false(any(grepl("--yes", line(), fixed = TRUE)))
  expect_true(any(grepl("--yes", line("--yes"), fixed = TRUE)))
})

test_that("unexported groups: a group with any stamped unit is not redone", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 1")
  u <- DBI::dbGetQuery(f$con, "SELECT ID, path, scaffold FROM annotate ORDER BY ID")
  assign_export_group(f$con, u, c("g1", "g1", "g2"))
  local_mocked_bindings(export_group_units = function(con, group) u)
  expect_equal(cli_unexported_groups(f$con), c("g1", "g2"))
  DBI::dbExecute(f$con, "UPDATE export SET export_time_stamp = 5 WHERE ID = 'A'")
  expect_equal(cli_unexported_groups(f$con), "g2")
})

test_that("group assignment respects existing members and never yields NA names", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 1")
  DBI::dbExecute(f$con, "UPDATE annotate SET partial = 'yes' WHERE ID IN ('B', 'C')")
  u <- DBI::dbGetQuery(f$con, "SELECT ID, path, scaffold FROM annotate ORDER BY ID")
  assign_export_group(f$con, u[1, ], "g")
  expect_error(cli_assign_group(f$con, u[2, ], "g", split = FALSE), "mix", class = "mitopilot_cli_error")
  expect_equal(cli_assign_group(f$con, u[2, ], "g", split = TRUE), "g-partial")
  expect_equal(cli_assign_group(f$con, u[3, ], "g", split = TRUE), "g-partial")
  orig <- fetch_export_data
  local_mocked_bindings(fetch_export_data = function(...) {
    x <- orig(...); x$completeness[x$ID == "C"] <- NA; x
  })
  expect_error(cli_assign_group(f$con, u[3, ], "h", split = TRUE), "C", class = "mitopilot_cli_error")
  expect_true(is.na(cli_assign_group(f$con, u[3, ], "h", split = TRUE, strict = FALSE)))
  expect_equal(auto_q(f, "SELECT COUNT(*) FROM export WHERE export_group = 'h'"), 0)
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

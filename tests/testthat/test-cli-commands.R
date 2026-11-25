spy <- function(name, env) {
  f <- function(...) assign("got", as.list(match.call())[-1], envir = env)
  formals(f) <- formals(get(name, envir = asNamespace("MitoPilot")))
  f
}

test_that("init --test assemblies calls the userAsmb test project with passthrough flags", {
  e <- environment()
  local_mocked_bindings(new_test_project_userAsmb = spy("new_test_project_userAsmb", e))
  expect_equal(cli_cmd_init(c("proj", "--test", "assemblies", "--n", "3", "--no-metadata", "--executor", "slurm")), 0L)
  expect_equal(e$got[c("path", "n", "fetch_metadata", "executor", "Rproj")],
               list(path = "proj", n = 3, fetch_metadata = FALSE, executor = "slurm", Rproj = FALSE))
})

test_that("init rejects an unknown test type and a flag the creator lacks", {
  expect_error(cli_cmd_init(c("p", "--test", "bogus")), class = "mitopilot_cli_error")
  expect_error(cli_cmd_init(c("p", "--test", "annotations", "--full-size")), class = "mitopilot_cli_error")
  expect_error(cli_cmd_init(c("p", "--test", "reads", "--n", "x")), class = "mitopilot_cli_error")
  expect_error(cli_cmd_init(c("p", "--settings", "a")), "opts", class = "mitopilot_cli_error")
})

test_that("init picks the creator and maps flags", {
  e <- environment()
  local_mocked_bindings(new_project = spy("new_project", e))
  expect_equal(cli_cmd_init(c("--mapping", "m.csv", "--data", "reads", "--min-depth", "5", "--genetic-code", "2", "--no-metadata")), 0L)
  got <- e$got
  expect_equal(got$mapping_fn, "m.csv"); expect_equal(got$data_path, "reads")
  expect_equal(got$min_depth, 5); expect_equal(got$genetic_code, 2)
  expect_false(got$fetch_geome); expect_false(got$Rproj)
  expect_error(cli_cmd_init(character()), class = "mitopilot_cli_error")
  expect_error(cli_cmd_init(c("--assemblies", "a", "--min-depth", "5")), class = "mitopilot_cli_error")
})

test_that("init --help prints usage and exits 0", {
  expect_message(st <- cli_cmd_init("--help"), "Usage: mitopilot init")
  expect_equal(st, 0L)
})

test_that("status --json emits only JSON and exits 2 when something is parked", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE, depth = c(80, 2))
  out <- capture.output(st <- suppressMessages(cli_cmd_status(c("--project", f$path, "--json"))))
  d <- jsonlite::fromJSON(paste(out, collapse = ""))
  expect_equal(st, 2L); expect_equal(d$reasons[d$ID == "B"], "low_coverage")
})

test_that("status table, filters, and exit code", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE, depth = c(80, 2))
  msg <- testthat::capture_messages(out <- capture.output(
    st <- cli_cmd_status(c("--project", f$path, "--parked"))))
  expect_equal(st, 2L)
  expect_true(any(grepl("B", out))); expect_false(any(grepl("^ *A ", out)))
  expect_match(msg[1], "2 units")
  capture.output(suppressMessages(st <- cli_cmd_status(c("--project", f$path, "--stage", "export"))))
  expect_equal(st, 2L)
})

test_that("status exits 0 when nothing is parked", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80)
  capture.output(suppressMessages(st <- cli_cmd_status(c("--project", f$path))))
  expect_equal(st, 0L)
})

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

q1 <- function(f, sql) DBI::dbGetQuery(f$con, sql)[[1]]

test_that("status rejects a bad --stage and stray arguments", {
  f <- local_project_fixture(ids = "A")
  expect_error(cli_cmd_status(c("--project", f$path, "--stage", "bogus")), class = "mitopilot_cli_error")
  expect_error(cli_cmd_status(c("extra", "--project", f$path)), class = "mitopilot_cli_error")
})

test_that("accept refuses a parked sample without --override and exits 2", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  expect_equal(suppressMessages(cli_cmd_accept(c("A", "--project", f$path))), 2L)
  expect_equal(q1(f, "SELECT assemble_lock FROM assemble"), 0)
})

test_that("accept --override needs --yes when scripted, then locks and logs", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("accept", "A", "--override", "--project", f$path))), 1L)
  expect_equal(suppressMessages(cli_main(c("accept", "A", "--override", "--yes", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT action FROM review_log"), "override")
})

test_that("accept of a passing sample locks it and logs nothing", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80)
  expect_equal(suppressMessages(cli_cmd_accept(c("A", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT assemble_lock FROM assemble"), 1)
  expect_equal(q1(f, "SELECT COUNT(*) FROM review_log"), 0)
})

test_that("accept --all-parked needs --override", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  expect_error(cli_cmd_accept(c("--all-parked", "--project", f$path)), "override", class = "mitopilot_cli_error")
  expect_equal(suppressMessages(cli_cmd_accept(c("--all-parked", "--override", "--yes", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT assemble_lock FROM assemble"), 1)
})

test_that("accept infers the annotate stage and refuses on warnings", {
  f <- local_project_fixture(ids = "A", annotated = TRUE, warnings = 2)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 0")
  expect_equal(suppressMessages(cli_cmd_accept(c("A", "--project", f$path))), 2L)
  expect_equal(suppressMessages(cli_cmd_accept(c("A:1:1", "--override", "--yes", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT annotate_lock FROM annotate"), 1)
})

test_that("accept asks for the unit spelling when a sample has several annotate units", {
  f <- local_project_fixture(ids = "A", annotated = TRUE,
    annotate_units = data.frame(ID = "A", path = 1, scaffold = 1:2))
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 0")
  expect_error(cli_cmd_accept(c("A", "--project", f$path)), "A:1:1", class = "mitopilot_cli_error")
})

test_that("hold sets state 0", {
  f <- local_project_fixture(ids = "A")
  expect_equal(suppressMessages(cli_cmd_hold(c("A", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT assemble_switch FROM assemble"), 0)
})

test_that("requeue needs --yes when scripted and changes nothing without it", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE)
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("requeue", "A", "--project", f$path))), 1L)
  expect_equal(q1(f, "SELECT assemble_switch FROM assemble"), 2)
  expect_equal(suppressMessages(cli_main(c("requeue", "A", "--yes", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT assemble_switch FROM assemble"), 1)
})

test_that("requeue unlocks a locked unit", {
  f <- local_project_fixture(ids = "A", annotated = TRUE)
  expect_equal(suppressMessages(cli_cmd_requeue(c("A:1:1", "--yes", "--project", f$path))), 0L)
  expect_equal(q1(f, "SELECT annotate_lock FROM annotate"), 0)
  expect_equal(q1(f, "SELECT annotate_switch FROM annotate"), 1)
})

test_that("run refuses while another run is active", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(run_active = function(...) "assemble_20261125_101500")
  expect_equal(suppressMessages(cli_main(c("run", "assemble", "--project", f$path))), 1L)
})

test_that("run --submit --dry-run prints a script and launches nothing", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  local_mocked_bindings(cli_self_command = function() "/x/mitopilot")
  local_mocked_bindings(run = function(...) stop("must not launch"), .package = "processx")
  out <- capture.output(st <- suppressMessages(cli_cmd_run(c("assemble", "--submit", "--dry-run", "--project", f$path))))
  expect_equal(st, 0L); expect_true(any(grepl("mitopilot run assemble", out, fixed = TRUE)))
  expect_false(any(grepl("--submit", grep("mitopilot run", out, value = TRUE), fixed = TRUE)))
  expect_true(any(grepl("--yes", grep("mitopilot run", out, value = TRUE), fixed = TRUE)))
})

test_that("run --dry-run prints the nextflow command", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(run_pipeline = function(...) stop("must not launch"))
  out <- capture.output(st <- suppressMessages(cli_cmd_run(c("assemble", "--dry-run", "--project", f$path))))
  expect_equal(st, 0L); expect_true(any(grepl("nextflow", out)))
})

test_that("run annotate fails while imported gene names are unconfirmed", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(annot_import_pending = function(...) TRUE)
  expect_error(cli_cmd_run(c("annotate", "--project", f$path)), "app", class = "mitopilot_cli_error")
})

test_that("run exits 1 on a failed run and 2 when units are parked", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(run_pipeline = function(...) list(status = 3L, base = NA_character_, log = "x.log", queued = 1L))
  expect_equal(suppressMessages(cli_cmd_run(c("assemble", "--project", f$path))), 1L)
  local_mocked_bindings(run_pipeline = function(...) {
    DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2, topology = 'circular'")
    list(status = 0L, base = NA_character_, log = NA_character_, queued = 1L)
  })
  DBI::dbExecute(f$con, "UPDATE assemblies SET depth = '1 1 1'")
  expect_equal(suppressMessages(cli_cmd_run(c("assemble", "--project", f$path))), 2L)
})

test_that("rules sets and validates the pass-rule switches", {
  f <- local_project_fixture(ids = "A")
  expect_equal(suppressMessages(cli_cmd_rules(c("--min-coverage", "25", "--accept-linear", "--project", f$path))), 0L)
  expect_equal(project_settings(f$con), list(min_mean_depth = 25, accept_linear = TRUE))
  expect_error(cli_cmd_rules(c("--min-coverage", "lots", "--project", f$path)), class = "mitopilot_cli_error")
  expect_error(cli_cmd_rules(c("--min-coverage", "-1", "--project", f$path)), class = "mitopilot_cli_error")
  suppressMessages(cli_cmd_rules(c("--no-accept-linear", "--project", f$path)))
  expect_false(project_settings(f$con)$accept_linear)
})

test_that("group rejects unlocked units and bad names", {
  f <- local_project_fixture(ids = "A", annotated = TRUE, warnings = 0)
  expect_error(cli_cmd_group(c("g 1", "--samples", "A", "--project", f$path)), class = "mitopilot_cli_error")
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 0")
  expect_error(cli_cmd_group(c("g1", "--samples", "A", "--project", f$path)), class = "mitopilot_cli_error")
})

test_that("group --split separates complete and partial genomes", {
  f <- local_project_fixture(ids = c("A", "B"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 1")
  DBI::dbExecute(f$con, "UPDATE annotate SET partial = 'yes' WHERE ID = 'B'")
  expect_error(cli_cmd_group(c("g", "--all-accepted", "--project", f$path)), "--split", class = "mitopilot_cli_error")
  expect_equal(suppressMessages(cli_cmd_group(c("g", "--all-accepted", "--split", "--project", f$path))), 0L)
  g <- DBI::dbGetQuery(f$con, "SELECT ID, export_group FROM export ORDER BY ID")
  expect_equal(g$export_group, c("g-complete", "g-partial"))
})

test_that("metadata fetch calls the chosen source", {
  f <- local_project_fixture(ids = "A")
  got <- character()
  local_mocked_bindings(
    fetch_geome = function(path = ".", ...) got <<- c(got, "geome"),
    fetch_gbif = function(path = ".", ...) got <<- c(got, "gbif"),
    fetch_ncbi = function(path = ".", ...) got <<- c(got, "ncbi"))
  suppressMessages(cli_cmd_metadata(c("fetch", "--source", "gbif", "--project", f$path)))
  expect_equal(got, "gbif")
  suppressMessages(cli_cmd_metadata(c("fetch", "--project", f$path)))
  expect_equal(got, c("gbif", "geome", "gbif", "ncbi"))
  expect_error(cli_cmd_metadata(c("fetch", "--source", "x", "--project", f$path)), class = "mitopilot_cli_error")
})

test_that("new commands print usage with --help", {
  for (cmd in c("run", "accept", "hold", "requeue", "group", "metadata", "rules")) {
    expect_message(st <- cli_main(c(cmd, "--help")), paste0("Usage: mitopilot ", cmd))
    expect_equal(st, 0L)
  }
})

test_that("accept names the refused sample and its reasons", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  msg <- testthat::capture_messages(cli_cmd_accept(c("A", "--project", f$path)))
  expect_true(any(grepl("A (low_coverage)", msg, fixed = TRUE)))
})

test_that("run --dry-run=true never launches or submits", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  local_mocked_bindings(run_pipeline = function(...) stop("must not launch"),
                        submit_pipeline = function(project, command, workflow, dry_run = FALSE) {
                          if (!dry_run) stop("must not submit")
                          list(script = "s.sh", lines = "x", submitted = FALSE, output = character())
                        })
  capture.output(st <- suppressMessages(cli_cmd_run(c("assemble", "--dry-run=true", "--submit=yes", "--project", f$path))))
  expect_equal(st, 0L)
})

test_that("run --submit hands submit_pipeline a clean command and reports failure", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  got <- NULL
  local_mocked_bindings(cli_self_command = function() "/x/mitopilot",
    submit_pipeline = function(project, command, workflow, dry_run = FALSE) {
      got <<- command
      list(script = "s.sh", lines = "", submitted = FALSE, output = "boom")
    })
  expect_equal(suppressMessages(cli_cmd_run(c("assemble", "--submit", "--project", f$path))), 1L)
  expect_false(grepl("--submit", got, fixed = TRUE))
  expect_true(grepl("--yes", got, fixed = TRUE))
  expect_true(grepl(paste0("--project '", normalizePath(f$path)), got, fixed = TRUE))
})

test_that("accept --override on mixed input writes nothing when the confirm fails", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE, depth = c(80, 2))
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("accept", "A", "B", "--override", "--project", f$path))), 1L)
  expect_equal(q1(f, "SELECT SUM(assemble_lock) FROM assemble"), 0)
})

test_that("a declined override is reported as declined", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  withr::local_options(MitoPilot.cli_interactive = TRUE, MitoPilot.cli_input = textConnection("n"))
  msg <- testthat::capture_messages(st <- cli_cmd_accept(c("A", "--override", "--project", f$path)))
  expect_equal(st, 2L)
  expect_true(any(grepl("declined", msg))); expect_false(any(grepl("Refused", msg)))
})

test_that("accept --reason needs --override and --stage is validated everywhere", {
  f <- local_project_fixture(ids = "A")
  expect_error(cli_cmd_accept(c("A", "--reason", "x", "--project", f$path)), "override", class = "mitopilot_cli_error")
  expect_error(cli_cmd_accept(c("--all-parked", "--override", "--stage", "bogus", "--project", f$path)), class = "mitopilot_cli_error")
  expect_error(cli_cmd_hold(c("A", "--stage", "bogus", "--project", f$path)), class = "mitopilot_cli_error")
  expect_error(cli_cmd_requeue(c("A", "--stage", "bogus", "--project", f$path)), class = "mitopilot_cli_error")
})

test_that("blocked units are not marked partial", {
  f <- local_project_fixture(ids = "A", assembled = TRUE, depth = 80, annotated = TRUE, topology = "linear")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbExecute(f$con, "UPDATE project_settings SET value = '1' WHERE key = 'accept_linear'")
  expect_equal(suppressMessages(cli_cmd_accept(c("A", "--stage", "assemble", "--project", f$path))), 2L)
  expect_equal(q1(f, "SELECT partial FROM annotate"), "no")
})

test_that("--all-parked override counts a sample once", {
  f <- local_project_fixture(ids = "A", assembled = TRUE, depth = 3,
    annotate_units = data.frame(ID = "A", path = 1, scaffold = 1:2))
  expect_equal(nrow(cli_parked_units(f$con, "assemble")), 1L)
})

test_that("init copies the saved submit template for --executor", {
  td <- withr::local_tempdir()
  prof <- file.path(td, "prof"); dir.create(prof)
  writeLines("#SBATCH x", file.path(prof, "slurm.submit.template"))
  proj <- file.path(td, "proj")
  local_mocked_bindings(
    new_test_project = function(path = ".", n = Inf, full_size = FALSE, executor = "local",
                                container = "", Rproj = FALSE, force = FALSE, fetch_metadata = TRUE) {
      dir.create(path)
    },
    mitopilot_config_dir = function() prof)
  expect_message(cli_cmd_init(c(proj, "--test", "reads", "--executor", "slurm")), "submit template")
  expect_true(file.exists(submit_template_path(proj)))
})

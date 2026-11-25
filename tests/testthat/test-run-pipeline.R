test_that("pipeline_queue matches the app's assemble gate", {
  f <- local_project_fixture(ids = c("A", "B", "C", "D"))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2 WHERE ID IN ('B','C','D')")
  DBI::dbExecute(f$con, "UPDATE assemble SET join_switch = 1 WHERE ID IN ('C','D')")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'D'")
  expect_setequal(pipeline_queue(f$con, "assemble"), c("A", "C"))
})

test_that("pipeline_queue user-assembly gate does not admit state 4", {
  f <- local_project_fixture(ids = c("A", "B"))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 4 WHERE ID = 'B'")
  expect_setequal(pipeline_queue(f$con, "assemble"), c("A", "B"))
  expect_equal(pipeline_queue(f$con, "assemble", user_asmb = TRUE), "A")
})

test_that("pipeline_queue annotate needs a locked sample and a ready unit", {
  f <- local_project_fixture(ids = c("A", "B"),
    annotate_units = data.frame(ID = c("A", "B"), path = 1, scaffold = 1))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'A'")
  expect_equal(pipeline_queue(f$con, "annotate"), "A")
})

test_that("annot_import_pending is TRUE only for unconfirmed gene names", {
  con <- local_test_db()
  expect_false(annot_import_pending(con))
  DBI::dbExecute(con, "CREATE TABLE gene_name_map (raw_name TEXT, confirmed INTEGER)")
  DBI::dbExecute(con, "INSERT INTO gene_name_map VALUES ('cox1', 0)")
  expect_true(annot_import_pending(con))
  DBI::dbExecute(con, "UPDATE gene_name_map SET confirmed = 1")
  expect_false(annot_import_pending(con))
})

test_that("nextflow_env matches the app's environment", {
  f <- local_project_fixture(ids = "A")
  pin <- nf_pin_version()
  expect_identical(nextflow_env(f$path), c(
    "current", if (!is.na(pin)) c(NXF_VER = pin), NXF_ANSI_SUMMARY = "TRUE", COLUMNS = "500",
    SGE = "/cm/shared/apps/uge/8.8.1/age", SGE_ARCH = "lx-amd64", SGE_CELL = "age",
    SGE_ROOT = "/cm/shared/apps/uge/8.8.1"))
})

test_that("run_pipeline does not start Nextflow on an empty queue", {
  f <- local_project_fixture(ids = "A")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2")
  local_mocked_bindings(run = function(...) stop("must not launch"), .package = "processx")
  r <- run_pipeline(f$path, "assemble", echo = FALSE)
  expect_equal(c(r$status, r$queued), c(0L, 0L))
})

test_that("run_pipeline runs Nextflow and writes a command-line run report", {
  f <- local_project_fixture(ids = "A")
  dir.create(file.path(f$path, ".nextflow", "cache"), recursive = TRUE)
  seen <- NULL
  local_mocked_bindings(run = function(command, args, wd, env, ...) {
    seen <<- list(command = command, args = args, wd = wd, env = env)
    log <- args[which(args == "-log") + 1]
    dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
    writeLines("Nov-25 10:00:00.000 [main] DEBUG nextflow.cli.Launcher - Execution complete -- Goodbye", log)
    list(status = 0L, stdout = "", stderr = "")
  }, .package = "processx")
  r <- run_pipeline(f$path, "assemble", resume = FALSE, echo = FALSE)
  expect_equal(c(r$status, r$queued), c(0L, 1L))
  expect_equal(seen$command, "nextflow")
  expect_equal(seen$wd, normalizePath(f$path))
  expect_identical(seen$env, nextflow_env(normalizePath(f$path)))
  expect_false("-resume" %in% seen$args)
  expect_true(all(c("-entry", "WF1") %in% seen$args))
  dcf <- read.dcf(run_report_paths(normalizePath(f$path), r$base)[["dcf"]])
  expect_equal(dcf[1, "launch"], c(launch = "cli"))
  expect_equal(dcf[1, "result"], c(result = "finished"))
})

test_that("run_pipeline uses the user-assembly entry and keeps -resume by default", {
  f <- local_project_fixture(ids = "A")
  dir.create(file.path(f$path, ".nextflow", "cache"), recursive = TRUE)
  args <- NULL
  local_mocked_bindings(run = function(command, args, ...) {
    args <<- args
    list(status = 1L, stdout = "", stderr = "")
  }, .package = "processx")
  r <- run_pipeline(f$path, "assemble", user_asmb = TRUE, echo = FALSE)
  expect_equal(r$status, 1L)
  expect_true(all(c("WF1_userAsmb", "-resume") %in% args))
})

test_that("run_pipeline refuses while a run is active", {
  f <- local_project_fixture(ids = "A")
  d <- run_dir(f$path, "nextflow")
  dir.create(d, recursive = TRUE)
  writeLines("running", file.path(d, paste0(run_basename("assemble"), ".nextflow.log")))
  local_mocked_bindings(run = function(...) stop("must not launch"), .package = "processx")
  expect_error(run_pipeline(f$path, "assemble", echo = FALSE), class = "mitopilot_cli_error")
})

test_that("submit_pipeline dry run writes a script holding the command and submits nothing", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  r <- submit_pipeline(f$path, "mitopilot auto --project /p --yes", "assemble", dry_run = TRUE)
  expect_false(r$submitted); expect_true(file.exists(r$script))
  expect_true(any(grepl("mitopilot auto --project /p --yes", r$lines, fixed = TRUE)))
  expect_true(any(grepl("^#SBATCH", r$lines)))
  expect_equal(dirname(r$script), run_dir(normalizePath(f$path), "jobs"))
})

test_that("submit_pipeline fills a saved template as is", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  writeLines(c("#!/bin/sh", "#SBATCH -J <<JOB_NAME>>", "<<NEXTFLOW_CMD>>"),
             submit_template_path(f$path))
  r <- submit_pipeline(f$path, "mitopilot auto --yes", "assemble", dry_run = TRUE)
  expect_equal(r$lines[3], "mitopilot auto --yes")
  expect_match(r$lines[2], "^#SBATCH -J assemble_")
})

test_that("submit_pipeline submits through run_submit", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  local_mocked_bindings(run_submit = function(script_path, executor, work_dir) {
    list(success = TRUE, command = "sbatch", output = paste("Submitted", executor))
  })
  r <- submit_pipeline(f$path, "mitopilot auto --yes", "assemble")
  expect_true(r$submitted)
  expect_equal(r$output, "Submitted slurm")
})

test_that("submit_pipeline refuses the local executor", {
  f <- local_project_fixture(ids = "A", executor = "local")
  expect_error(submit_pipeline(f$path, "mitopilot auto", "assemble", dry_run = TRUE), class = "mitopilot_cli_error")
})

local_saved_profile <- function(name, scheduler = "slurm", env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  suppressMessages(generate_config(name, scheduler = scheduler, container_engine = "singularity",
                                   queue = "short", profile_dir = d))
  writeLines(cluster_build_template(scheduler, "short"), file.path(d, paste0(name, ".submit.template")))
  testthat::local_mocked_bindings(mitopilot_config_dir = function() d, .env = env)
  d
}

ok_tools <- function(...) c(sbatch = TRUE, nextflow = TRUE, java = TRUE, singularity = TRUE)
ok_nf <- function(...) list(ok = TRUE, detail = "nextflow x")
ok_cfg <- function(...) list(ok = TRUE, detail = "")

local_stub_bin <- function(tools = c("nextflow", "java"), env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  for (t in tools) {
    writeLines(c("#!/bin/sh", "echo stub version 1"), file.path(d, t))
    Sys.chmod(file.path(d, t), "0755")
  }
  withr::local_path(d, action = "prefix", .local_envir = env)
  d
}

fake_submit <- function(script_path, executor, work_dir) {
  s <- readLines(script_path)
  log <- sub("^#SBATCH -o ", "", grep("^#SBATCH -o ", s, value = TRUE))
  system2("sh", shQuote(script_path), stdout = log, stderr = log)
  list(success = TRUE, command = "sbatch", output = "Submitted batch job 1")
}

test_that("cluster_check stops at the first failing step and reports the tool's error", {
  local_saved_profile("mycl", scheduler = "slurm")
  local_mocked_bindings(cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg,
                        run_submit = function(...) list(success = FALSE, command = "sbatch",
                          output = "sbatch: error: invalid partition specified: bigmem"))
  r <- cluster_check("mycl")
  expect_equal(r$step, c("tools", "config", "submit"))
  expect_equal(r$ok, c(TRUE, TRUE, FALSE))
  expect_match(r$detail[3], "invalid partition", fixed = TRUE)
})

test_that("cluster_check stops at tools and names what is missing", {
  local_saved_profile("mycl")
  local_mocked_bindings(cluster_tools_found = function(...) c(sbatch = FALSE, nextflow = TRUE, java = FALSE))
  r <- cluster_check("mycl")
  expect_equal(r$step, "tools")
  expect_match(r$detail, "sbatch, java", fixed = TRUE)
})

test_that("cluster_check passes submit when the job log reports success, then runs the pipeline step", {
  local_saved_profile("mycl")
  local_stub_bin()
  local_mocked_bindings(cluster_head_tools = function() c("nextflow", "java"))
  local_mocked_bindings(cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg, run_submit = fake_submit,
                        cluster_pipeline_runs = function(...) list(ok = TRUE, detail = "loads"))
  r <- cluster_check("mycl", wait = 20)
  expect_equal(r$step, c("tools", "config", "submit", "pipeline"))
  expect_true(all(r$ok))
})

test_that("cluster_check times out when no job log appears", {
  local_saved_profile("mycl")
  local_mocked_bindings(cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg,
                        run_submit = function(...) list(success = TRUE, command = "sbatch", output = "Submitted batch job 9"))
  cancelled <- NULL
  local_mocked_bindings(cluster_cancel = function(executor, id) cancelled <<- c(executor, id))
  r <- cluster_check("mycl", wait = 0)
  expect_false(r$ok[3])
  expect_match(r$detail[3], "timed out after 0 s", fixed = TRUE)
  expect_match(r$detail[3], "Submitted batch job 9", fixed = TRUE)
  expect_equal(cancelled, c("slurm", "9"))
})

test_that("a template missing a required token fails the submit stage without submitting", {
  d <- local_saved_profile("mycl")
  writeLines(c("#SBATCH -J x", "nextflow run foo"), file.path(d, "mycl.submit.template"))
  called <- FALSE
  local_mocked_bindings(cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg,
                        run_submit = function(...) { called <<- TRUE; list(success = TRUE, output = "1") })
  r <- cluster_check("mycl")
  expect_false(r$ok[3])
  expect_match(r$detail[3], "<<NEXTFLOW_CMD>> and <<LOG_FILE>>", fixed = TRUE)
  expect_false(called)
})

test_that("the check job reports missing tools from inside the job", {
  local_saved_profile("mycl")
  local_stub_bin()
  local_mocked_bindings(cluster_head_tools = function() c("nextflow", "java", "no_such_tool_xyz"),
                        cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg, run_submit = fake_submit)
  r <- cluster_check("mycl", wait = 20)
  expect_false(r$ok[3])
  expect_match(r$detail[3], "missing: no_such_tool_xyz", fixed = TRUE)
})

test_that("the check job fails when nextflow does not start", {
  local_saved_profile("mycl")
  d <- local_stub_bin()
  writeLines(c("#!/bin/sh", "echo boom >&2", "exit 3"), file.path(d, "nextflow"))
  local_mocked_bindings(cluster_head_tools = function() c("nextflow", "java"),
                        cluster_tools_found = ok_tools, cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg, run_submit = fake_submit)
  r <- cluster_check("mycl", wait = 20)
  expect_false(r$ok[3])
  expect_match(r$detail[3], "nextflow -version failed: boom", fixed = TRUE)
})

test_that("a stage that throws becomes a failed row", {
  local_saved_profile("mycl")
  local_mocked_bindings(cluster_tools_found = function(...) stop("cannot launch"))
  r <- cluster_check("mycl")
  expect_equal(r$step, "tools"); expect_false(r$ok)
  expect_match(r$detail, "cannot launch")
})

test_that("a skipped submit stage is reported as skipped, not passed", {
  f <- local_project_fixture(executor = "local")
  local_mocked_bindings(cluster_tools_found = function(...) c(nextflow = TRUE), cluster_nextflow_version = ok_nf,
                        cluster_config_parses = ok_cfg,
                        cluster_pipeline_runs = function(...) list(ok = TRUE, detail = "loads"))
  msgs <- paste(testthat::capture_messages(st <- cli_cmd_cluster(c("check", "--project", f$path))), collapse = "")
  expect_equal(st, 0L)
  expect_match(msgs, "submit +skipped")
  expect_match(msgs, "3 passed, 1 skipped")
  expect_false(grepl("All checks passed", msgs))
})

test_that("check rejects both a NAME and --project", {
  local_saved_profile("mycl")
  expect_error(cli_cmd_cluster(c("check", "mycl", "--project", ".")), "not both", class = "mitopilot_cli_error")
})

test_that("cluster_check without a name uses the project's own config and template", {
  f <- local_project_fixture(executor = "slurm")
  seen <- NULL
  local_mocked_bindings(cluster_tools_found = function(ctx) { seen <<- ctx; ok_tools() },
                        cluster_nextflow_version = ok_nf, cluster_config_parses = ok_cfg,
                        run_submit = function(...) list(success = FALSE, command = "sbatch", output = "no"))
  writeLines("#SBATCH -J <<JOB_NAME>>\n<<NEXTFLOW_CMD>>", file.path(f$path, ".mitopilot_submit.template"))
  r <- cluster_check(project = f$path)
  expect_equal(seen$executor, "slurm")
  expect_match(seen$template, ".mitopilot_submit.template", fixed = TRUE)
  expect_equal(r$ok, c(TRUE, TRUE, FALSE))
  expect_false(any(grepl("mitopilot_check", list.files(f$path, all.files = TRUE))))
})

test_that("setup writes a profile and a tokenized submit template with the head-job answers", {
  d <- withr::local_tempdir()
  withr::local_options(MitoPilot.cli_interactive = TRUE,
    MitoPilot.cli_input = textConnection(c("mycl", "1", "short", "", "", "1", "/scratch/cache",
                                           "24:00:00", "8G", "module load nextflow", "", "n")))
  local_mocked_bindings(mitopilot_config_dir = function() d)
  expect_equal(suppressMessages(cli_cmd_cluster("setup")), 0L)
  tmpl <- readLines(file.path(d, "mycl.submit.template"))
  expect_true(all(c("#SBATCH -p short", "#SBATCH -t 24:00:00", "#SBATCH --mem=8G", "module load nextflow") %in% tmpl))
  expect_true(any(grepl("<<NEXTFLOW_CMD>>", tmpl, fixed = TRUE)))
  expect_true(any(grepl("<<LOG_FILE>>", tmpl, fixed = TRUE)))
  expect_true("mycl" %in% list_configs(d)$name)
  expect_match(paste(readLines(file.path(d, "config.mycl")), collapse = "\n"), "cacheDir = '/scratch/cache'", fixed = TRUE)
})

test_that("template builder sets time, memory, account, and setup lines per scheduler", {
  p <- cluster_build_template("pbs", "q", account = "acct", time = "10:00:00", mem = "8G")
  expect_true(all(c("#PBS -l select=1:ncpus=1:mem=8gb", "#PBS -l walltime=10:00:00", "#PBS -A acct") %in% p))
  s <- cluster_build_template("sge", penv = "smp", time = "10:00:00", mem = "8G")
  expect_true(all(c("#$ -pe smp 1", "#$ -l h_rt=10:00:00", "#$ -l h_vmem=8G",
                    'export NXF_OPTS="-Xms256m -Xmx4096m"') %in% s))
  expect_lt(match('export NXF_OPTS="-Xms256m -Xmx4096m"', s), match("<<NEXTFLOW_CMD>>", s))
  s2 <- cluster_build_template("sge", mem = "lots")
  expect_false(any(grepl("h_vmem|NXF_OPTS", s2)))
  l <- cluster_build_template("lsf", time = "10:00:00", mem = "8G")
  expect_true(all(c("#BSUB -M 8000", "#BSUB -W 10:00") %in% l))
  n <- cluster_build_template("slurm", native_prefix = "/opt/mp", setup = "module load x")
  expect_true(any(grepl("activate.sh", n, fixed = TRUE)))
  expect_lt(match("module load x", n), match("<<NEXTFLOW_CMD>>", n))
  expect_false(any(grepl("EDIT: load", cluster_build_template("slurm", setup = "m"))))
})

test_that("a saved template is accepted unchanged by the project script builder", {
  d <- local_saved_profile("mycl")
  f <- local_project_fixture(executor = "slurm")
  expect_true(cluster_copy_submit_template("mycl", f$path))
  s <- build_submit_script(f$path, "slurm", "short", "mitopilot run assemble", "job1", "/x/job1.log")
  expect_true(all(c("#SBATCH -J job1", "#SBATCH -o /x/job1.log", "mitopilot run assemble") %in% s))
  expect_false(any(grepl("<<", s, fixed = TRUE)))
})

test_that("cluster_copy_submit_template copies only when a saved template exists and the project has none", {
  d <- withr::local_tempdir(); local_mocked_bindings(mitopilot_config_dir = function() d)
  proj <- withr::local_tempdir()
  expect_false(cluster_copy_submit_template("mycl", proj))
  writeLines("X <<NEXTFLOW_CMD>>", file.path(d, "mycl.submit.template"))
  expect_true(cluster_copy_submit_template("mycl", proj))
  expect_equal(readLines(file.path(proj, ".mitopilot_submit.template")), "X <<NEXTFLOW_CMD>>")
  writeLines("mine", file.path(proj, ".mitopilot_submit.template"))
  expect_false(cluster_copy_submit_template("mycl", proj))
  expect_equal(readLines(file.path(proj, ".mitopilot_submit.template")), "mine")
})

test_that("cluster setup and edit refuse to run without a terminal", {
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("cluster", "setup"))), 1L)
  local_saved_profile("mycl")
  expect_equal(suppressMessages(cli_main(c("cluster", "edit", "mycl"))), 1L)
})

test_that("cluster list and show print the saved files", {
  local_saved_profile("mycl")
  expect_output(expect_equal(cli_cmd_cluster("list"), 0L), "mycl")
  out <- capture.output(cli_cmd_cluster(c("show", "mycl")))
  expect_true(any(grepl("config.mycl", out, fixed = TRUE)))
  expect_true(any(grepl("<<NEXTFLOW_CMD>>", out, fixed = TRUE)))
  expect_error(cli_cmd_cluster("show"), class = "mitopilot_cli_error")
  expect_error(cli_cmd_cluster("bogus"), class = "mitopilot_cli_error")
})

test_that("cluster apply updates a project's config and template", {
  local_saved_profile("mycl")
  f <- local_project_fixture(executor = "local")
  cfg <- file.path(f$path, ".config")
  writeLines(fill_config(readLines(cfg), list(CONTAINER_ID = "img", RAW_DIR = "r", ASMB_DIR = "NA",
                                              MIN_DEPTH = "5", NCBI_API_KEY = "", CONTAINER_ENGINE = "docker { enabled = true }")), cfg)
  expect_equal(suppressMessages(cli_cmd_cluster(c("apply", "mycl", "--project", f$path, "--yes"))), 0L)
  expect_equal(read_config_executor(file.path(f$path, ".config"))$executor, "slurm")
  expect_true(file.exists(file.path(f$path, ".mitopilot_submit.template")))
})

test_that("cluster apply backs up an existing project template and keeps it when the profile has none", {
  d <- local_saved_profile("mycl")
  f <- local_project_fixture(executor = "local")
  cfg <- file.path(f$path, ".config")
  writeLines(fill_config(readLines(cfg), list(CONTAINER_ID = "img", RAW_DIR = "r", ASMB_DIR = "NA",
                                              MIN_DEPTH = "5", NCBI_API_KEY = "", CONTAINER_ENGINE = "docker { enabled = true }")), cfg)
  tp <- file.path(f$path, ".mitopilot_submit.template")
  writeLines("old <<NEXTFLOW_CMD>>", tp)
  suppressMessages(cli_cmd_cluster(c("apply", "mycl", "--project", f$path, "--yes")))
  bak <- list.files(f$path, "^\\.mitopilot_submit\\.template\\.bak\\.", all.files = TRUE, full.names = TRUE)
  expect_length(bak, 1)
  expect_equal(readLines(bak), "old <<NEXTFLOW_CMD>>")
  expect_false(identical(readLines(tp), "old <<NEXTFLOW_CMD>>"))
  unlink(file.path(d, "mycl.submit.template"))
  writeLines("mine", tp)
  suppressMessages(cli_cmd_cluster(c("apply", "mycl", "--project", f$path, "--yes")))
  expect_equal(readLines(tp), "mine")
})

test_that("cluster check prints a hint and exits 1 on failure", {
  local_saved_profile("mycl")
  local_mocked_bindings(cluster_check = function(...) data.frame(step = c("tools", "config"), ok = c(TRUE, FALSE),
                                                                 detail = c("x", "parse error")))
  msgs <- testthat::capture_messages(st <- cli_cmd_cluster(c("check", "mycl")))
  expect_equal(st, 1L)
  expect_match(paste(msgs, collapse = ""), "parse error")
  expect_match(paste(msgs, collapse = ""), "mitopilot cluster edit mycl", fixed = TRUE)
})

# MitoPilot CLI Mode, Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A `mitopilot` shell command that takes a project from `init` to GenBank export unattended for clean samples, parks the rest with reasons, and can be set up and smoke-tested on a cluster.

**Architecture:** A shell wrapper hands off to `MitoPilot:::cli_main()`, which parses flags and dispatches to one `cli_cmd_*` function per subcommand. The commands and the Shiny observers both call a new Shiny-free shared core (state, lock, export group, run launch, pass rules) extracted from app server code. Stage 2 (`opts`) and stage 3 (`review`, focused app launch) get their own plans.

**Tech Stack:** R (>= 4.4), DBI/RSQLite, dplyr, processx, cli, jsonlite, testthat, chromote (app regression), Nextflow 24.10.x to 25.10.x.

**Spec:** `tools/cli_mode_design.md`

## Global Constraints

- No new package dependencies (`cli`, `processx`, `jsonlite` already imported).
- No non-ASCII characters in code. Minimal comments.
- Nextflow starts with the project directory as working directory; version checks and pin (`check_nextflow_version`, `nf_pin_version`, `native_nf_pin`) unchanged.
- App behaviour for users who never use the command line must not change. Every task that edits an `app_*.R` file ends with the browser click-through in Task 12's protocol for the actions it touched.
- Exit codes: 0 done and nothing parked; 1 error; 2 finished but units parked or flagged.
- Any command that makes finished units rerun, or overrides a pass rule, needs `--yes` when stdin is not a terminal.
- Pass-rule thresholds: mean coverage default `10`; outlier defaults `start_aa = 10`, `stop_aa = 10`, `ident_pct = 60`.
- A parked or flagged unit never advances without an explicit override, and every override is logged.
- Work on branch `cli-mode`. Commit per task, brief messages, no attribution lines, never push.
- New NEWS.md entries go under `# MitoPilot (development version)`.

## Review Focus

1. Command run outside a project directory, or with `--project` pointing at a non-project: clear one-line error, exit 1, no stack trace, nothing created.
2. A pipeline run already active on the project when `run` or `auto` starts: refuse with a message naming the active run; never two pipelines on one project.
3. Pass-rule inputs that are missing rather than bad (NULL/NA coverage, NA warnings, sample never run): unit is parked with a reason, never passed and never an R error.
4. Old projects whose database lacks the new tables and columns: `open_project()` migrates them before any command reads them.
5. Unknown subcommand, unknown flag, or flag missing its value: usage text for that command, exit 1, nothing written.

---

## File Structure

| File | Responsibility |
|---|---|
| `inst/cli/mitopilot` | shell wrapper template (`@RSCRIPT@` placeholder) |
| `inst/cli/check.nf` | one-task pipeline for `cluster check` |
| `R/cli_main.R` | `cli_main`, `cli_parse_args`, usage text, exit codes |
| `R/cli_prompt.R` | `cli_is_interactive`, `cli_readline`, `cli_ask`, `cli_confirm` |
| `R/cli_install.R` | `install_cli` (exported) |
| `R/cli_project.R` | `open_project`, `project_settings`, `set_project_setting`, `run_active` |
| `R/project_actions.R` | shared core: state, locks, export group |
| `R/run_pipeline.R` | shared core: queue, env, launch, submit |
| `R/pass_rules.R` | pass rules, `auto_accept`, `log_review`, `project_status` |
| `R/cli_commands.R` | `cli_cmd_init/status/run/accept/hold/requeue/group/metadata/export/auto` |
| `R/cli_cluster.R` | `cli_cmd_cluster`, `cluster_check`, saved submit templates |
| `R/init_db.R`, `R/init_db_userAsmb.R`, `R/backwards_compatibility.R` | new tables and columns |
| `R/app_assemble.R`, `R/app_assemble_utils.R`, `R/app_annotate.R`, `R/app_export.R`, `R/app_run_pipline.R`, `R/app_run_pipline_userAsmb.R` (and userAsmb twins of the above where they duplicate the logic) | observers call the shared core |

Unit type used throughout: `units` is a data.frame with columns `ID`, `path`, `scaffold`. Assembly-stage functions take `ids` (character).

---

### Task 1: Flag reader, dispatcher, prompts

**Files:** Create `R/cli_main.R`, `R/cli_prompt.R`, `tests/testthat/test-cli-main.R`

**Interfaces:**
- Produces:
  - `cli_parse_args(args, flags = character()) -> list(pos = character(), opts = named list)`. `--key value` and `--key=value` give strings; names in `flags` give `TRUE`; `--no-<flag>` gives `FALSE`; dashes in names become underscores; tokens after `--` are positional. Unknown `--no-x` where `x` is not in `flags`, or a value flag with no value, calls `cli_fail()`.
  - `cli_check_opts(opts, known)`: `cli_fail()` naming any option not in `known`. Every `cli_cmd_*` calls it before doing anything else.
  - `cli_fail(msg, status = 1L)`: signals a condition of class `mitopilot_cli_error` carrying `status`.
  - `cli_main(args = commandArgs(trailingOnly = TRUE)) -> integer` exit status. Dispatches `args[1]` to `cli_cmd_<name>(args[-1])`; each command returns its integer status. Catches `mitopilot_cli_error` (message to stderr, its status) and any other error (message to stderr, 1). `--help`, no args, or unknown command print usage; unknown command returns 1.
  - `cli_is_interactive() -> logical` (`isatty(stdin())`, overridable by `options(MitoPilot.cli_interactive=)`).
  - `cli_readline(prompt) -> character(1)`: reads one line from `getOption("MitoPilot.cli_input", file("stdin"))`.
  - `cli_ask(prompt, choices) -> character(1)`: numbered menu over a named character vector (names are returned keys); re-asks on bad input.
  - `cli_confirm(prompt, yes = FALSE) -> logical`: returns `TRUE` when `yes`; when not interactive and not `yes`, `cli_fail("... needs --yes")`.

- [ ] **Step 1: Write failing tests**

```r
test_that("cli_parse_args splits positionals, values, and flags", {
  p <- cli_parse_args(c("assemble", "--no-resume", "--project", "/x", "--min-depth=5", "--submit"),
                      flags = c("resume", "submit"))
  expect_equal(p$pos, "assemble")
  expect_identical(p$opts, list(resume = FALSE, project = "/x", min_depth = "5", submit = TRUE))
})
test_that("value flag without a value fails with status 1", {
  e <- tryCatch(cli_parse_args("--project"), mitopilot_cli_error = identity)
  expect_s3_class(e, "mitopilot_cli_error"); expect_equal(e$status, 1L)
})
test_that("cli_main returns 1 and writes nothing for unknown command", {
  expect_equal(suppressMessages(cli_main("nonsense")), 1L)
})
test_that("cli_check_opts rejects an unknown option", {
  expect_error(cli_check_opts(list(project = "/x", bogus = TRUE), c("project", "json")),
               class = "mitopilot_cli_error")
})
test_that("cli_ask re-asks until a valid choice", {
  con <- textConnection(c("9", "2")); withr::local_options(MitoPilot.cli_input = con)
  expect_equal(suppressMessages(cli_ask("Pick", c(a = "A", b = "B"))), "b")
})
test_that("cli_confirm fails without --yes when not interactive", {
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_error(cli_confirm("Rerun 3 samples?"), class = "mitopilot_cli_error")
  expect_true(cli_confirm("Rerun 3 samples?", yes = TRUE))
})
```

- [ ] **Step 2:** Run `Rscript -e 'devtools::test(filter = "cli-main")'`; expect failures (functions missing).
- [ ] **Step 3:** Implement the functions above. Prompts and menus go to stderr via `message()` so stdout stays clean for `--json`.
- [ ] **Step 4:** Rerun; expect PASS.
- [ ] **Step 5:** Commit `cli: flag reader, dispatcher, prompts`.

### Task 2: Shell wrapper and install paths

**Files:** Create `inst/cli/mitopilot`, `R/cli_install.R`, `tests/testthat/test-cli-install.R`. Modify `inst/native/install_mitopilot_native.sh` (after the main env is built), `docker/Dockerfile` (after `install_local`), `NAMESPACE` via roxygen.

**Interfaces:**
- Produces: `install_cli(dir = "~/.local/bin", rscript = file.path(R.home("bin"), "Rscript"), overwrite = TRUE) -> invisible(path)` (exported). Writes the wrapper with `@RSCRIPT@` replaced, mode 0755, creates `dir`, and messages a PATH hint when `dir` is not on `PATH`.
- Wrapper body (exact):

```sh
#!/bin/sh
exec "@RSCRIPT@" --vanilla -e 'quit(status = MitoPilot:::cli_main())' "$@"
```

`cli_main()` default args already use `commandArgs(trailingOnly = TRUE)`.

- [ ] **Step 1: Failing test**

```r
test_that("install_cli writes an executable wrapper bound to an Rscript", {
  d <- withr::local_tempdir()
  p <- suppressMessages(install_cli(dir = d, rscript = "/opt/R/bin/Rscript"))
  expect_true(file.exists(p)); expect_equal(basename(p), "mitopilot")
  expect_true(any(grepl('exec "/opt/R/bin/Rscript"', readLines(p), fixed = TRUE)))
  expect_equal(unname(file.access(p, 1)), 0)
})
```

- [ ] **Step 2:** Run, expect fail. **Step 3:** Implement `install_cli`.
- [ ] **Step 4:** Native installer: write the same wrapper to `$main/bin/mitopilot` with `$main/bin/Rscript`. Dockerfile: `RUN Rscript -e 'MitoPilot::install_cli(dir = "/usr/local/bin")'`.
- [ ] **Step 5:** Verify by hand: `Rscript -e 'devtools::install(quick = TRUE)'; Rscript -e 'MitoPilot::install_cli(dir = "/tmp/mpbin")'; /tmp/mpbin/mitopilot --help; echo $?` prints usage and `0`; `/tmp/mpbin/mitopilot nonsense; echo $?` prints `1`.
- [ ] **Step 6:** Commit `cli: wrapper script and install_cli`.

### Task 3: Project handle, new tables, migration

**Files:** Create `R/cli_project.R`, `tests/testthat/test-cli-project.R`. Modify `R/init_db.R`, `R/init_db_userAsmb.R`, `R/backwards_compatibility.R`.

**Schema (both `new_db*` and a migration step that adds whatever is missing):**

```sql
CREATE TABLE project_settings (key TEXT PRIMARY KEY, value TEXT);
-- seeded: ('min_mean_depth','10'), ('accept_linear','0')
CREATE TABLE review_log (
  ID TEXT NOT NULL, path INTEGER, scaffold INTEGER,
  stage TEXT NOT NULL,          -- 'assemble' | 'annotate' | 'export'
  action TEXT NOT NULL,         -- 'auto' | 'override'
  rules TEXT,                   -- comma-separated failing reason codes, '' for auto
  reason TEXT, time_stamp INTEGER);
ALTER TABLE export_opts ADD COLUMN start_aa INTEGER;      -- default row: 10
ALTER TABLE export_opts ADD COLUMN stop_aa INTEGER;       -- 10
ALTER TABLE export_opts ADD COLUMN ident_pct REAL;        -- 60
ALTER TABLE export_opts ADD COLUMN aa_alignments INTEGER; -- 1
ALTER TABLE export_opts ADD COLUMN gene_export INTEGER;   -- 0
```

**Interfaces:**
- Produces:
  - `open_project(path = ".", migrate = TRUE) -> list(path, db, con, dir_out, user_asmb)`. `path` normalized; `cli_fail()` if `.sqlite` or `.config` is missing; sets `options(MitoPilot.db = db)`; runs `backwards_compatibility()` the way `app_server` does when `migrate`; `dir_out` derived exactly as `R/app_server.R:98-103`; `user_asmb` is `TRUE` when `.config` `asmbDir` is not `"NA"` (same test as `run_app`). Caller disconnects `con`.
  - `project_settings(con) -> named list` with `min_mean_depth` numeric and `accept_linear` logical.
  - `set_project_setting(con, key, value)`; unknown key is `cli_fail()`.
  - `run_active(path) -> character(1)` run base name of a live run or `NA`. A run is live when the newest `.runs/nextflow/*.nextflow.log` has no report under `.runs/reports` and `run_progress_text()` does not report it stale.

- [ ] **Step 1: Failing tests** (build the DB with `new_db()` on a temp mapping fixture, as `tests/testthat/test-backwards-compatibility.R` does):

```r
test_that("new projects carry settings, review_log, and export option defaults", {
  con <- local_test_db()   # helper in this file wrapping new_db on a tempdir
  expect_equal(project_settings(con), list(min_mean_depth = 10, accept_linear = FALSE))
  expect_true(DBI::dbExistsTable(con, "review_log"))
  eo <- DBI::dbGetQuery(con, "SELECT * FROM export_opts WHERE export_opts = 'default'")
  expect_equal(c(eo$start_aa, eo$stop_aa, eo$ident_pct, eo$aa_alignments, eo$gene_export), c(10, 10, 60, 1, 0))
})
test_that("migration adds the new pieces to an old database", {
  con <- local_test_db()
  DBI::dbExecute(con, "DROP TABLE project_settings"); DBI::dbExecute(con, "DROP TABLE review_log")
  DBI::dbExecute(con, "ALTER TABLE export_opts DROP COLUMN ident_pct")
  migrate_cli_schema(con)
  expect_equal(project_settings(con)$min_mean_depth, 10)
  expect_true("ident_pct" %in% DBI::dbListFields(con, "export_opts"))
})
test_that("open_project fails cleanly outside a project", {
  expect_error(open_project(withr::local_tempdir()), class = "mitopilot_cli_error")
})
test_that("set_project_setting rejects unknown keys", {
  con <- local_test_db()
  expect_error(set_project_setting(con, "bogus", 1), class = "mitopilot_cli_error")
})
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement; `migrate_cli_schema(con)` is the idempotent migration, called from `backwards_compatibility()` and listed wherever `schema_gaps` enumerates expected tables.
- [ ] **Step 4:** Run this file plus `devtools::test(filter = "backwards")`; PASS.
- [ ] **Step 5:** Commit `cli: project handle, settings and review log tables`.

### Task 4: Shared core for state, locks, export group; app calls it

**Files:** Create `R/project_actions.R`, `tests/testthat/test-project-actions.R`. Modify `R/app_assemble.R:598-624`, `R/app_assemble_utils.R` (`assemble_lock_begin`, `apply_assemble_lock`), `R/app_annotate.R:741-870` (`update_state`, `write_lock`, `toggle_lock`), `R/app_export.R` (`assign_export_group`, `make_group` name check), and the userAsmb twins that repeat the same writes.

**Interfaces:**
- Produces:
  - `set_unit_state(con, stage = c("assemble", "annotate"), units, state) -> data.frame` of rows written. `state` in 0:4. For `assemble`, `units` needs only `ID`.
  - `lock_assemblies(con, dir_out, ids, lock = TRUE) -> list(written = data.frame(ID, assemble_lock, ...), stale = data.frame)`. Locking: IDs returned by `stale_assemble_dirs(con, dir_out, ids = ids, pending_only = FALSE)` are not locked and are returned in `stale`; `join_switch` set `NA`; `time_stamp` set. Unlocking: no stale check, `join_switch` untouched.
  - `lock_annotations(con, units, lock = TRUE) -> list(written = data.frame, multi_path = character())`. Locking is all-or-nothing: if, after applying, any sample would have more than one distinct `path` locked, write nothing and return those IDs in `multi_path`.
  - `valid_group_name(name) -> logical` (`^[A-Za-z0-9._-]+$`).
  - `assign_export_group(con, units, group) -> data.frame` upserted into `export`; `group` is one name, a vector aligned to `units`, or `NA` to clear.
- App observers keep their dialogs, toasts, and `rv` updates; only the database write and the block logic move. `rv$data` is updated from the returned `written` rows.

- [ ] **Step 1: Failing tests** (database fixture from Task 3; insert rows directly):

```r
test_that("lock_assemblies holds back samples whose output is not on disk", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = "A")   # B has no output dir
  r <- lock_assemblies(f$con, f$dir_out, c("A", "B"))
  expect_equal(r$written$ID, "A"); expect_equal(r$stale$ID, "B")
  got <- DBI::dbGetQuery(f$con, "SELECT ID, assemble_lock, join_switch FROM assemble ORDER BY ID")
  expect_equal(got$assemble_lock, c(1, 0)); expect_true(is.na(got$join_switch[1]))
})
test_that("lock_annotations refuses two paths of one sample and writes nothing", {
  f <- local_project_fixture(ids = "A", annotate_units = data.frame(ID = "A", path = 1:2, scaffold = 1))
  r <- lock_annotations(f$con, data.frame(ID = "A", path = 1:2, scaffold = 1))
  expect_equal(r$multi_path, "A"); expect_equal(nrow(r$written), 0)
  expect_equal(sum(DBI::dbGetQuery(f$con, "SELECT annotate_lock FROM annotate")$annotate_lock, na.rm = TRUE), 0)
})
test_that("lock_annotations counts already locked paths", {
  f <- local_project_fixture(ids = "A", annotate_units = data.frame(ID = "A", path = 1:2, scaffold = 1))
  lock_annotations(f$con, data.frame(ID = "A", path = 1, scaffold = 1))
  expect_equal(lock_annotations(f$con, data.frame(ID = "A", path = 2, scaffold = 1))$multi_path, "A")
})
test_that("set_unit_state writes only the named units", {
  f <- local_project_fixture(ids = c("A", "B"))
  set_unit_state(f$con, "assemble", data.frame(ID = "A"), 0)
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_switch FROM assemble ORDER BY ID")[[1]], c(0, 1))
})
test_that("assign_export_group upserts and clears", {
  f <- local_project_fixture(ids = "A", annotate_units = data.frame(ID = "A", path = 1, scaffold = 1))
  u <- data.frame(ID = "A", path = 1, scaffold = 1)
  assign_export_group(f$con, u, "g1")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT export_group FROM export")[[1]], "g1")
  assign_export_group(f$con, u, NA)
  expect_true(is.na(DBI::dbGetQuery(f$con, "SELECT export_group FROM export")[[1]]))
  expect_false(valid_group_name("bad name")); expect_true(valid_group_name("fish_2026-A.1"))
})
```

Put `local_project_fixture()` in `tests/testthat/helper-cli.R`; later tasks reuse it.

- [ ] **Step 2:** Run, fail. **Step 3:** Implement `R/project_actions.R`.
- [ ] **Step 4:** Before touching app code, record the baseline: run the Task 12 click-through for state, lock, unlock, and export group on the unmodified app and save the database row dumps.
- [ ] **Step 5:** Change the observers to call the new functions. Run `devtools::test()` (whole suite) and expect PASS, then rerun the click-through and diff the row dumps against the baseline; expect no differences apart from `time_stamp`.
- [ ] **Step 6:** Commit `core: shared state, lock, and export group actions`.

### Task 5: Shared run launcher; app calls it

**Files:** Create `R/run_pipeline.R`, `tests/testthat/test-run-pipeline.R`. Modify `R/app_run_pipline.R:140-163, 325-364`, `R/app_run_pipline_userAsmb.R` (same spots), `R/app_annot_import.R` (split the check out of `annot_import_blocks_run`).

**Interfaces:**
- Consumes: `nextflow_cmd()`, `run_basename()`, `run_dir()`, `run_report()`, `write_run_report()`, `read_config_executor()`, `native_env()`, `nf_pin_version()`, `native_nf_pin()`, `build_submit_script()`, `run_submit()`, `run_progress_text()`.
- Produces:
  - `pipeline_queue(con, workflow = c("assemble", "annotate")) -> character` IDs, with the exact filters now in `app_run_pipline.R:140-163`.
  - `annot_import_pending(con) -> logical`: the non-Shiny half of `annot_import_blocks_run`.
  - `nextflow_env(project) -> character` named env vector exactly as built at `app_run_pipline.R:342-362` (including `"current"` first).
  - `run_pipeline(project, workflow, resume = TRUE, user_asmb = FALSE, echo = TRUE) -> list(status = integer, base = character, log = character, queued = integer)`. Fails with `cli_fail()` if `run_active(project)` is not `NA`. Returns `status = 0L, queued = 0L` without starting Nextflow when the queue is empty. Otherwise runs `processx::run("nextflow", args, wd = project, env = nextflow_env(project), error_on_status = FALSE)`, streaming output when `echo`, then writes the run report (`how = "cli"`).
  - `submit_pipeline(project, command, workflow, dry_run = FALSE) -> list(script = character path, lines = character, submitted = logical, output = character)`. Builds the job script with `build_submit_script(full_nf_cmd = command, ...)` where `command` is the `mitopilot ...` command line to run inside the job; writes it under `.runs/jobs/<base>.sh`; submits with `run_submit()` unless `dry_run`. Executor `local` is `cli_fail("no scheduler configured; see 'mitopilot cluster setup'")`.
- The app keeps its own `processx::process$new()` (it needs the live handle for the progress modal) but takes `args`' environment from `nextflow_env()` and the queue from `pipeline_queue()`.

- [ ] **Step 1: Failing tests**

```r
test_that("pipeline_queue matches the app's assemble gate", {
  f <- local_project_fixture(ids = c("A", "B", "C", "D"))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2 WHERE ID IN ('B','C','D')")
  DBI::dbExecute(f$con, "UPDATE assemble SET join_switch = 1 WHERE ID IN ('C','D')")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'D'")
  expect_setequal(pipeline_queue(f$con, "assemble"), c("A", "C"))
})
test_that("pipeline_queue annotate needs a locked sample and a ready unit", {
  f <- local_project_fixture(ids = c("A", "B"),
    annotate_units = data.frame(ID = c("A", "B"), path = 1, scaffold = 1))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'A'")
  expect_equal(pipeline_queue(f$con, "annotate"), "A")
})
test_that("run_pipeline does not start Nextflow on an empty queue", {
  f <- local_project_fixture(ids = "A")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 2")
  local_mocked_bindings(run = function(...) stop("must not launch"), .package = "processx")
  r <- run_pipeline(f$path, "assemble", echo = FALSE)
  expect_equal(c(r$status, r$queued), c(0L, 0L))
})
test_that("submit_pipeline dry run writes a script holding the command and submits nothing", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  r <- submit_pipeline(f$path, "mitopilot auto --project /p --yes", "assemble", dry_run = TRUE)
  expect_false(r$submitted); expect_true(file.exists(r$script))
  expect_true(any(grepl("mitopilot auto --project /p --yes", r$lines, fixed = TRUE)))
  expect_true(any(grepl("^#SBATCH", r$lines)))
})
test_that("submit_pipeline refuses the local executor", {
  f <- local_project_fixture(ids = "A", executor = "local")
  expect_error(submit_pipeline(f$path, "mitopilot auto", "assemble", dry_run = TRUE), class = "mitopilot_cli_error")
})
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement.
- [ ] **Step 4:** Baseline click-through of "Update" on the unmodified app for both workflows on the reads test project (Task 12 protocol), then switch the two `pipeline_server` files to the shared pieces and repeat; a run must start, finish, and write its report identically. Whole suite PASS.
- [ ] **Step 5:** Commit `core: shared pipeline queue, environment, launch, and submit`.

### Task 6: Pass rules, auto-accept, status

**Files:** Create `R/pass_rules.R`, `tests/testthat/test-pass-rules.R`.

**Interfaces:**
- Consumes: `lock_assemblies`, `lock_annotations`, `project_settings`.
- Produces:
  - `assembly_pass(con, ids = NULL, settings = project_settings(con)) -> data.frame(ID, pass, reasons)`; `reasons` is a comma-separated string of codes, `""` when passing. Codes and tests, all evaluated (no short-circuit):
    - `not_run`: `assemble_switch` in 0, 1, 4 or NA
    - `failed`: `assemble_switch == 3`
    - `multiple_paths`: more than one non-ignored `path` in `assemblies`
    - `multiple_scaffolds`: more than one non-ignored `scaffold`
    - `linear`: `assemble.topology` is not `"circular"`, unless `settings$accept_linear`
    - `low_coverage`: mean per-base depth of the non-ignored unit below `settings$min_mean_depth`, or depth missing. Parse `assemblies.depth` the way the caller of `score_assembly_paths()` does at `R/app_assemble_coverage_details.R:118`; lift that parsing into a small shared function rather than copying it.
    - `poor_reference`: `assemble.poor_blast_ref` is set to its "yes" value
  - `annotation_pass(con, units = NULL) -> data.frame(ID, path, scaffold, pass, reasons)`. Codes: `not_run` (`annotate_switch` in 0, 1, NA), `failed` (3), `warnings` (`annotate.warnings` is NA or > 0), `problematic` (`annotate.problematic == "yes"`). Only units whose sample has `assemble_lock == 1` are considered.
  - `export_flags(project, group) -> list(flagged = data.frame(ID, gene, why), result = <flag_PCG_outliers value>, report = character path)`. Thresholds from the `default` row of `export_opts`. Always renders the HTML alignment report via `make_PCG_alignments()` into `<dir_out>/export/<group>/` and returns its path. `flagged` has one row per sample and gene with `why` in `start_offset`, `stop_offset`, `low_identity`, `internal_stop`.
  - `log_review(con, units, stage, action = c("auto", "override"), rules = "", reason = NA)`.
  - `auto_accept(con, dir_out, stage = c("assemble", "annotate")) -> list(accepted = data.frame, parked = data.frame)`. Locks units that pass and are not yet locked, logs them with `action = "auto"`. With `accept_linear`, a passing linear unit also gets `annotate.partial = "yes"` once its annotate row exists. Units held back by the lock functions (`stale`, `multi_path`) are returned as parked with reason `output_missing` or `multiple_paths`.
  - `override_accept(con, dir_out, stage, units, reason = NA) -> list(accepted, blocked)`. Locks regardless of pass rules, logs `action = "override"` with the failing codes. The lock functions' own blocks still apply and come back in `blocked`.
  - `project_status(con) -> data.frame(ID, path, scaffold, stage, state, locked, accepted_by, reasons, export_group)`; `stage` is the furthest stage reached; `accepted_by` is `"auto"`, `"override"`, `"manual"`, or `NA` (latest `review_log` row for that unit and stage; a lock with no log row is `"manual"`).
  - `warning_tally(con) -> data.frame(warning, n)` over `annotations.warnings`.

- [ ] **Step 1: Failing tests.** One test per reason code plus:

```r
test_that("clean circular single-path sample passes and is locked as auto", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80)
  expect_true(assembly_pass(f$con)$pass)
  r <- auto_accept(f$con, f$dir_out, "assemble")
  expect_equal(r$accepted$ID, "A")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT action FROM review_log")[[1]], "auto")
  expect_equal(project_status(f$con)$accepted_by, "auto")
})
test_that("missing depth parks instead of erroring", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = NA)
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
})
test_that("several failing rules are all reported", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3, topology = "linear")
  expect_setequal(strsplit(assembly_pass(f$con)$reasons, ",")[[1]], c("linear", "low_coverage"))
})
test_that("accept_linear passes a clean linear sample", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80, topology = "linear")
  set_project_setting(f$con, "accept_linear", 1)
  expect_true(assembly_pass(f$con)$pass)
})
test_that("any curation warning parks; NA warnings parks", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE, warnings = c(0, 1, NA))
  expect_equal(annotation_pass(f$con)$pass, c(TRUE, FALSE, FALSE))
  expect_equal(annotation_pass(f$con)$reasons, c("", "warnings", "warnings"))
})
test_that("auto_accept never locks a parked unit; override does and logs the rules", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  expect_equal(nrow(auto_accept(f$con, f$dir_out, "assemble")$accepted), 0)
  override_accept(f$con, f$dir_out, "assemble", data.frame(ID = "A"), reason = "re-sequencing not possible")
  lg <- DBI::dbGetQuery(f$con, "SELECT action, rules, reason FROM review_log")
  expect_equal(c(lg$action, lg$rules), c("override", "low_coverage"))
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_lock FROM assemble")[[1]], 1)
})
```

Extend `local_project_fixture()` with the `assembled`, `depth`, `topology`, `annotated`, `warnings` arguments used here.

- [ ] **Step 2:** Run, fail. **Step 3:** Implement. **Step 4:** PASS.
- [ ] **Step 5:** Commit `core: pass rules, auto-accept, overrides, status`.

### Task 7: `init`, `status`

**Files:** Create `R/cli_commands.R`, `tests/testthat/test-cli-commands.R`.

**Interfaces:**
- Produces `cli_cmd_init(args)`, `cli_cmd_status(args)`, each returning an exit status.
- `init [DIR]` flags: `--test reads|assemblies|annotations`, `--mapping`, `--data`, `--assemblies DIR`, `--annotations DIR`, `--n`, `--full-size`, `--no-raw-data`, `--executor`, `--container`, `--force`, `--no-metadata`, `--min-depth`, `--genetic-code`. Mapping:
  - `--test reads` -> `new_test_project`; `assemblies` -> `new_test_project_userAsmb`; `annotations` -> `new_test_project_userAnnot`. Any other value is `cli_fail()` listing the three.
  - no `--test`: `--annotations` -> `new_project_userAnnot`; `--assemblies` -> `new_project_userAsmb`; else `new_project` (needs `--mapping` and `--data`).
  - `--no-metadata` -> `fetch_metadata = FALSE`; `Rproj = FALSE` always. Flags a chosen creator does not accept are `cli_fail()` (check against `names(formals())`), not silently dropped.
  - `--settings` is reserved: `cli_fail("--settings arrives with 'mitopilot opts'")`.
- `status` flags: `--project`, `--json`, `--stage`, `--parked`. Prints `project_status()` as a table (to stdout) followed by one summary line and, when any annotation warnings exist, the `warning_tally()`. `--json` prints `jsonlite::toJSON(..., dataframe = "rows", na = "null")` and nothing else on stdout. Exit 2 when any unit is parked, else 0.

- [ ] **Step 1: Failing tests** (mock the creators with `local_mocked_bindings` and capture their arguments):

```r
test_that("init --test assemblies calls the userAsmb test project with passthrough flags", {
  got <- NULL
  local_mocked_bindings(new_test_project_userAsmb = function(...) got <<- list(...))
  expect_equal(cli_cmd_init(c("proj", "--test", "assemblies", "--n", "3", "--no-metadata", "--executor", "slurm")), 0L)
  expect_equal(got[c("path", "n", "fetch_metadata", "executor", "Rproj")],
               list(path = "proj", n = 3, fetch_metadata = FALSE, executor = "slurm", Rproj = FALSE))
})
test_that("init rejects an unknown test type and a flag the creator lacks", {
  expect_error(cli_cmd_init(c("p", "--test", "bogus")), class = "mitopilot_cli_error")
  expect_error(cli_cmd_init(c("p", "--test", "annotations", "--full-size")), class = "mitopilot_cli_error")
})
test_that("status --json emits only JSON and exits 2 when something is parked", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE, depth = c(80, 2))
  out <- capture.output(st <- suppressMessages(cli_cmd_status(c("--project", f$path, "--json"))))
  d <- jsonlite::fromJSON(paste(out, collapse = ""))
  expect_equal(st, 2L); expect_equal(d$reasons[d$ID == "B"], "low_coverage")
})
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement. **Step 4:** PASS.
- [ ] **Step 5:** Commit `cli: init and status`.

### Task 8: `run`, `accept`, `hold`, `requeue`, `group`, `metadata`

**Files:** Modify `R/cli_commands.R`, `tests/testthat/test-cli-commands.R`.

**Interfaces:**
- `run assemble|annotate` flags: `--project`, `--no-resume`, `--submit`, `--dry-run`, `--auto-accept`, `--yes`.
  - `--submit`: `submit_pipeline()` with the same command line minus `--submit`, plus `--yes` and an absolute `--project`. `--dry-run` prints the script and the Nextflow command and exits 0 without running or submitting.
  - No `--submit` and `read_config_executor()$executor` is not `local` and no scheduler job variable is set (`SLURM_JOB_ID`, `JOB_ID`, `PBS_JOBID`, `LSB_JOBID`): warn that the head process will run on this node, then continue.
  - `annotate` with `annot_import_pending(con)`: `cli_fail()` telling the user to confirm gene names in the app.
  - After a run, `--auto-accept` calls `auto_accept()` for that stage. Exit: Nextflow status non-zero -> 1; else 2 if anything is parked at that stage; else 0.
- `accept SAMPLE... [--stage assemble|annotate] [--override] [--reason TEXT] [--all-parked] [--yes]`. Without `--override`, a unit failing its pass rules is refused with its reasons (exit 2) and passing units are locked (logged nowhere, so they show as `manual`). With `--override`, `cli_confirm()` then `override_accept()`. `--all-parked` requires `--override`. A sample name with several annotate units accepts the single non-ignored unit, or fails listing the units and the `ID:path:scaffold` form to use.
- `hold SAMPLE...` -> `set_unit_state(..., 0)`; `requeue SAMPLE...` -> unlock if locked, then state 1, behind `cli_confirm("N finished units will rerun")`.
- `group NAME --samples a,b | --all-accepted`: validates with `valid_group_name()`; units must have `annotate_lock == 1`; refuses a mix of complete and partial genomes unless `--split`, which makes `NAME-complete` and `NAME-partial` (the app's rule at `app_export.R` `make_group`).
- `metadata fetch [--source geome|gbif|ncbi]`: calls `fetch_geome` / `fetch_gbif` / `fetch_ncbi` for the project; default all three.
- `rules [--min-coverage N] [--accept-linear | --no-accept-linear]`: with no flags prints the current pass-rule settings; with flags calls `set_project_setting()` (`min_mean_depth`, `accept_linear`). A non-numeric or negative `--min-coverage` is `cli_fail()`. Add `mitopilot rules` to the spec's command list.

- [ ] **Step 1: Failing tests**

```r
test_that("accept refuses a parked sample without --override and exits 2", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  expect_equal(suppressMessages(cli_cmd_accept(c("A", "--project", f$path))), 2L)
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_lock FROM assemble")[[1]], 0)
})
test_that("accept --override needs --yes when scripted, then locks and logs", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("accept", "A", "--override", "--project", f$path))), 1L)
  expect_equal(suppressMessages(cli_main(c("accept", "A", "--override", "--yes", "--project", f$path))), 0L)
  expect_equal(DBI::dbGetQuery(f$con, "SELECT action FROM review_log")[[1]], "override")
})
test_that("run refuses while another run is active", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(run_active = function(...) "assemble_20261125_101500")
  expect_equal(suppressMessages(cli_main(c("run", "assemble", "--project", f$path))), 1L)
})
test_that("run --submit --dry-run prints a script and launches nothing", {
  f <- local_project_fixture(ids = "A", executor = "slurm")
  local_mocked_bindings(run = function(...) stop("must not launch"), .package = "processx")
  out <- capture.output(st <- suppressMessages(cli_cmd_run(c("assemble", "--submit", "--dry-run", "--project", f$path))))
  expect_equal(st, 0L); expect_true(any(grepl("mitopilot run assemble", out, fixed = TRUE)))
  expect_false(any(grepl("--submit", grep("mitopilot run", out, value = TRUE), fixed = TRUE)))
})
test_that("rules sets and validates the pass-rule switches", {
  f <- local_project_fixture(ids = "A")
  expect_equal(suppressMessages(cli_cmd_rules(c("--min-coverage", "25", "--accept-linear", "--project", f$path))), 0L)
  expect_equal(project_settings(f$con), list(min_mean_depth = 25, accept_linear = TRUE))
  expect_error(cli_cmd_rules(c("--min-coverage", "lots", "--project", f$path)), class = "mitopilot_cli_error")
})
test_that("group rejects unlocked units and bad names", {
  f <- local_project_fixture(ids = "A", annotated = TRUE, warnings = 0)
  expect_error(cli_cmd_group(c("g 1", "--samples", "A", "--project", f$path)), class = "mitopilot_cli_error")
  expect_error(cli_cmd_group(c("g1", "--samples", "A", "--project", f$path)), class = "mitopilot_cli_error")
})
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement. **Step 4:** PASS.
- [ ] **Step 5:** Commit `cli: run, accept, hold, requeue, group, metadata`.

### Task 9: `export` with flagged-sample handling; app export window uses saved options

**Files:** Modify `R/cli_commands.R`, `R/app_export.R` (export window inputs and `run_export`), `R/app_export_utils.R` (`get_export_opts` / `set_export_opts`), `tests/testthat/test-cli-export.R`.

**Interfaces:**
- Consumes: `export_flags()`, `export_files()`, `export_group_units()`, `log_review()`, `specimen_conflicts()`.
- `export GROUP [--flagged ask|stop|drop|keep] [--yes]`. Default `ask` when interactive, `stop` otherwise; `ask` when not interactive is `cli_fail()`.
  1. Blocks, each exit 1 with the app's wording: unknown or empty group; a sample with more than one path in the group; invalid header templates.
  2. Warnings printed, not blocking: specimen metadata conflicts, NMNH problems, fragmented samples (the app's "export anyway" confirmations). In a script they require `--yes`.
  3. `export_flags()`; print the report path every time.
  4. Nothing flagged -> export -> exit 0.
  5. Flagged: list `ID  gene  why`, then by choice:
     - `stop`: export nothing, exit 2.
     - `drop`: `export_files(IDs = <unflagged IDs>)`, exit 2.
     - `keep`: `log_review(stage = "export", action = "override", rules = <why codes>)` per flagged unit, export all, exit 0.
     - `ask`: menu "drop flagged / keep all / cancel" (the "review in app" entry is added in stage 3), then as above; cancel is exit 2.
  `export_files()` is always called with `review = FALSE`, `out_dir = dir_out`, and thresholds and switches from the `default` row of `export_opts`.
- App: the export window's start, stop, identity, alignment, and gene-export inputs are initialised from `export_opts` and written back when the user exports, through `get_export_opts` / `set_export_opts` extended with the five columns.

- [ ] **Step 1: Failing tests** (mock `export_flags` and `export_files`):

```r
flagged_two <- function(...) list(
  flagged = data.frame(ID = c("B", "B"), gene = c("ND5", "COX1"), why = c("low_identity", "internal_stop")),
  result = NULL, report = "/p/out/export/g1/AA_alignments_g1.html")
test_that("scripted export stops by default when samples are flagged", {
  f <- local_export_fixture(group = "g1", ids = c("A", "B"))
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  local_mocked_bindings(export_flags = flagged_two, export_files = function(...) stop("must not export"))
  msgs <- capture.output(st <- cli_cmd_export(c("g1", "--project", f$path)), type = "message")
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
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement the command. **Step 4:** PASS.
- [ ] **Step 5:** Baseline click-through of the app export on the docs export project, wire the export window to the saved options, repeat, and diff the written files (`.fasta`, `.tbl`, summary CSV byte-identical) and the `export` rows.
- [ ] **Step 6:** Commit `cli: export with flagged-sample choices; export options saved`.

### Task 10: `auto`

**Files:** Modify `R/cli_commands.R`, `tests/testthat/test-cli-auto.R`.

**Interfaces:**
- `auto [--group NAME] [--flagged stop|drop|keep] [--submit] [--dry-run] [--no-resume] [--yes]`. `--group` defaults to the project directory name. `--submit` and `--dry-run` behave as in `run`, submitting one job that runs `mitopilot auto ... --yes`.
- Chain, each step skipped when it has nothing to do:
  1. `run_pipeline(assemble)`; non-zero status stops with exit 1.
  2. `auto_accept(assemble)`.
  3. `run_pipeline(annotate)`; skipped with a message when `annot_import_pending()`.
  4. `auto_accept(annotate)`.
  5. Assign every accepted, ungrouped unit to `--group` (split complete and partial as in `group --split`).
  6. `export` for each resulting group with the given `--flagged` (default `stop`).
  7. Print the `status` summary. Exit 2 if anything is parked or was flagged, else 0.
- Rerunning `auto` after clearing parked units picks up where it left off (the chain is driven by database state, not by remembered progress).

- [ ] **Step 1: Failing tests** (mock `run_pipeline`, `cli_cmd_export`):

```r
test_that("auto accepts clean units, parks the rest, exports, and exits 2", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = c("A", "B"), assembled = TRUE,
                             depth = c(80, 2), annotated = TRUE, warnings = c(0, 0))
  calls <- character(); exported <- NULL
  local_mocked_bindings(
    run_pipeline = function(project, workflow, ...) { calls <<- c(calls, workflow); list(status = 0L, queued = 1L) },
    cli_cmd_export = function(args) { exported <<- args[1]; 0L })
  expect_equal(suppressMessages(cli_cmd_auto(c("--group", "g1", "--project", f$path, "--yes"))), 2L)
  expect_equal(calls, c("assemble", "annotate")); expect_equal(exported, "g1")
  st <- project_status(f$con)
  expect_equal(st$export_group[st$ID == "A"], "g1"); expect_true(is.na(st$export_group[st$ID == "B"]))
})
test_that("auto stops with exit 1 when the assemble run fails", {
  f <- local_project_fixture(ids = "A")
  local_mocked_bindings(run_pipeline = function(...) list(status = 1L, queued = 1L))
  expect_equal(suppressMessages(cli_cmd_auto(c("--project", f$path, "--yes"))), 1L)
})
```

- [ ] **Step 2:** Run, fail. **Step 3:** Implement. **Step 4:** PASS.
- [ ] **Step 5:** Commit `cli: auto`.

### Task 11: `cluster`

**Files:** Create `R/cli_cluster.R`, `inst/cli/check.nf`, `tests/testthat/test-cli-cluster.R`. Modify `R/generate_config.R` only if `resolve_config()` needs to return the template path.

**Interfaces:**
- Consumes: `generate_config()`, `list_configs()`, `mitopilot_config_dir()`, `resolve_config()`, `migrate_config()`, `submission_script()`, `submit_template_path()`, `tokenize_submit_script()`, `run_submit()`.
- Saved submit template lives beside the saved profile: `<mitopilot_config_dir()>/<name>.submit.template`, tokenized like the per-project one. `cli_cmd_init` copies it to the project's `.mitopilot_submit.template` when the chosen `--executor` has one (add this to Task 7's code here, with a test).
- Produces:
  - `cli_cmd_cluster(args)` with subcommands:
    - `setup [--name N]`: interactive only. Asks, in order: profile name; scheduler (`slurm`, `sge`, `pbs`, `lsf`); queue; account or project code (blank to skip); extra scheduler options (blank to skip); for `sge`, parallel environment (default `mthread`); container engine (`singularity`, `apptainer`, `docker`, `none`); container cache folder, or native install folder when `none`; head-job time limit (default `168:00:00`); head-job memory (default `16G`); setup lines to run before the pipeline (blank-terminated, e.g. `module load ...`). Calls `generate_config()` with the answers, builds the submit template from `submission_script()` with the head-job answers and setup lines substituted, saves it, then offers to run `cluster check`.
    - `list`, `show NAME` (prints both files), `edit NAME [--submit]` (opens the config, or the submit template with `--submit`, via `utils::file.edit` honouring `$EDITOR`; not interactive is `cli_fail()`), `apply NAME --project DIR` (`migrate_config()` plus copying the submit template, behind `cli_confirm`).
    - `check [NAME] [--project DIR]`.
  - `cluster_check(name = NULL, project = NULL, wait = 300) -> data.frame(step, ok, detail)`, stopping at the first failure. Steps:
    1. `tools`: scheduler submit command for the executor, `nextflow` (and `check_nextflow_version()`), `java`, and the container engine or native `activate.sh` are found.
    2. `config`: the profile parses (`cp` to a temp dir as `nextflow.config`, `nextflow config .`).
    3. `submit`: submit a job script built from the saved template whose command is `hostname`; poll for its log file up to `wait` seconds; `detail` is the scheduler's output on failure.
    4. `pipeline`: `nextflow run inst/cli/check.nf -c <profile>` in a temp dir; `check.nf` has one process, using the profile's container or native activation, whose script is `Rscript -e 'library(MitoPilot)'`.
    `detail` carries the tool's own error text. A fixed hint per step is printed by the command (for example step 3: "check queue name, account, and time limit in 'mitopilot cluster edit NAME --submit'").
  Exit 0 when all four pass, else 1.

- [ ] **Step 1: Failing tests** (no real scheduler; mock `Sys.which`, `processx::run`, `run_submit`):

```r
test_that("cluster_check stops at the first failing step and reports the tool's error", {
  local_saved_profile("mycl", scheduler = "slurm")     # helper: generate_config into a temp profile dir
  local_mocked_bindings(cluster_tools_found = function(...) c(sbatch = TRUE, nextflow = TRUE, java = TRUE, singularity = TRUE),
                        cluster_config_parses = function(...) list(ok = TRUE, detail = ""),
                        run_submit = function(...) list(success = FALSE, command = "sbatch",
                                                         output = "sbatch: error: invalid partition specified: bigmem"))
  r <- cluster_check("mycl")
  expect_equal(r$step, c("tools", "config", "submit")); expect_equal(r$ok, c(TRUE, TRUE, FALSE))
  expect_match(r$detail[3], "invalid partition", fixed = TRUE)
})
test_that("setup writes a profile and a tokenized submit template with the head-job answers", {
  d <- withr::local_tempdir(); withr::local_options(MitoPilot.cli_interactive = TRUE,
    MitoPilot.cli_input = textConnection(c("mycl", "1", "short", "", "", "1", "/scratch/cache", "24:00:00", "8G", "module load nextflow", "", "n")))
  local_mocked_bindings(mitopilot_config_dir = function() d)
  expect_equal(suppressMessages(cli_cmd_cluster("setup")), 0L)
  tmpl <- readLines(file.path(d, "mycl.submit.template"))
  expect_true(all(c("#SBATCH -p short", "#SBATCH -t 24:00:00", "#SBATCH --mem=8G", "module load nextflow") %in% tmpl))
  expect_true(any(grepl("<<NEXTFLOW_CMD>>", tmpl, fixed = TRUE)))
  expect_true("mycl" %in% list_configs(d)$name)
})
test_that("init copies a saved submit template into the new project", {
  d <- withr::local_tempdir(); writeLines("X <<NEXTFLOW_CMD>>", file.path(d, "mycl.submit.template"))
  local_mocked_bindings(mitopilot_config_dir = function() d,
                        new_test_project = function(path, ...) dir.create(path))
  p <- file.path(withr::local_tempdir(), "proj")
  cli_cmd_init(c(p, "--test", "reads", "--executor", "mycl"))
  expect_equal(readLines(file.path(p, ".mitopilot_submit.template")), "X <<NEXTFLOW_CMD>>")
})
test_that("cluster setup refuses to run without a terminal", {
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_equal(suppressMessages(cli_main(c("cluster", "setup"))), 1L)
})
```

Adjust the `list_configs()` assertion to its real return shape.

- [ ] **Step 2:** Run, fail. **Step 3:** Implement. **Step 4:** PASS.
- [ ] **Step 5:** `mitopilot run assemble --submit --dry-run` for each of slurm, sge, pbs, lsf on a scratch project: the script must equal what the app's "Write Submission Script" produces for the same project except for the command line.
- [ ] **Step 6:** Commit `cli: cluster setup, check, and saved submit templates`.

### Task 12: End-to-end runs and app regression

**Files:** Create `dev/cli_mode/e2e_cli.sh`, `dev/cli_mode/app_regression.R`, `dev/cli_mode/README.md` (how to run both; `dev/` is gitignored for outputs, scripts are committed if the directory's siblings are).

**App click-through protocol (used by Tasks 4, 5, 9 and rerun in full here).** Follow the existing harness notes (`dev/reorient_e2e` as the example; run `backwards_compatibility` on the scratch copy first; expect the container popup; use the docs export project as the export baseline). Drive the real app in a browser with chromote (`tests/testthat/helper-chromote.R`), on scratch copies of a reads project and a user-assembly project. For each action: click it in the UI, then dump the affected rows with `export_db_to_csv()` and compare with the baseline taken on `devel` at the commit the branch started from.

Actions, both app flavours unless noted:
1. Assemble: set state for one row and for a multi-row selection.
2. Assemble: lock one, lock several, lock with a missing output directory (expect the "Assembly output not found" alert and the rest locked), unlock with confirmation.
3. Annotate: set state; lock one path; try to lock two paths of one sample (expect the "Only one assembly path" alert, nothing written); unlock.
4. Export: create a group, split a mixed complete and partial selection, remove from group.
5. Update: launch assemble and annotate from the app on the small reads test project; run finishes; run report appears; progress modal behaves.
6. Export: export a group with and without flagged samples; export window values persist after reopening the window.

- [ ] **Step 1:** Write `app_regression.R` to run the list above and print a pass or fail line per action with the row diff on failure.
- [ ] **Step 2:** Run it against the baseline commit and against the branch head; expect every action identical apart from time stamps. Fix any difference in the task that introduced it.
- [ ] **Step 3:** Write `e2e_cli.sh`: for each of `reads`, `assemblies`, `annotations`: `mitopilot init --test <type> --n 3 <scratch>`, `mitopilot auto --project <scratch> --yes; echo $?`, `mitopilot status --project <scratch> --json`. Expected: exit 0 or 2 (never 1); every unit either exported or listed with a reason; exported `.fasta` and `.tbl` exist for the group; a second `mitopilot auto` on the same project launches no run and returns the same status.
- [ ] **Step 4:** Run it locally (docker executor) and record the three status tables in the README.
- [ ] **Step 5:** Review Focus checks by hand: `mitopilot status` in `/tmp` (exit 1, one line); `mitopilot run assemble` while another `mitopilot run` is active in a second terminal (refused); `mitopilot status` on a project copy made before this branch (migrates, then prints); `mitopilot run --bogus` (usage, exit 1).
- [ ] **Step 6:** If a cluster is available, `mitopilot cluster setup`, `cluster check`, and `mitopilot auto --submit` on the reads test project there. If not, record that it was not done.
- [ ] **Step 7:** Commit `cli: end-to-end and app regression harness`.

### Task 13: Documentation

**Files:** Create `vignettes/Command-Line.Rmd`. Modify `_pkgdown.yml` (articles list), `README.Rmd` (one short section and link), `vignettes/Custom-HPC.Rmd` (point to `mitopilot cluster`), `vignettes/Native-Install.Rmd` and `vignettes/Installation.Rmd` (how the command gets on PATH), `NEWS.md`.

- [ ] **Step 1:** Write the vignette: install the command (both install types), the three test projects, `auto`, reading `status` and parked reasons, overriding, export choices, cluster setup and check, exit codes for scripts. State that settings editing and the review walk-through are coming, and that parked samples are handled in the app meanwhile.
- [ ] **Step 2:** `Rscript -e 'devtools::document(); devtools::check(vignettes = FALSE)'`; expect no new errors, warnings, or notes.
- [ ] **Step 3:** Commit `docs: command-line mode`.

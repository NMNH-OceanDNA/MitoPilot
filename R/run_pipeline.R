#' IDs a run would update: samples in Assemble, one per unit in Annotate
#' @noRd
pipeline_queue <- function(con, workflow = c("assemble", "annotate"), user_asmb = FALSE) {
  workflow <- match.arg(tolower(workflow[1]), c("assemble", "annotate"))
  if (workflow == "assemble") {
    # join_switch = 1 is a join-only redo that WF1 admits on its own
    q <- dplyr::tbl(con, "assemble")
    q <- if (user_asmb) {
      dplyr::filter(q, assemble_switch == 1 | (join_switch == 1 & assemble_lock == 0))
    } else {
      dplyr::filter(q, assemble_switch %in% c(1, 4) | (join_switch == 1 & assemble_lock == 0))
    }
    return(dplyr::pull(q, ID))
  }
  dplyr::left_join(dplyr::tbl(con, "assemble"), dplyr::tbl(con, "annotate"), by = "ID") |>
    dplyr::filter(assemble_lock == 1 & annotate_switch == 1) |>
    dplyr::pull(ID)
}

#' TRUE when imported gene names still need confirming
#' @noRd
annot_import_pending <- function(con) {
  is_user_annot(con) && !user_annot_confirmed(con)
}

#' Environment for a Nextflow run of a project
#' @noRd
nextflow_env <- function(project) {
  nat <- read_config_executor(file.path(project, ".config"))$native_activate
  nat_env <- if (!is.null(nat)) native_env(nat) else character()
  nxf_pin <- if (length(nat_env)) native_nf_pin(nat_env) else nf_pin_version()
  c(
    "current",
    nat_env,
    if (!is.na(nxf_pin)) c(NXF_VER = nxf_pin),
    NXF_ANSI_SUMMARY = TRUE,
    # Keep Nextflow's ANSI log from truncating process names
    COLUMNS = 500,
    SGE = "/cm/shared/apps/uge/8.8.1/age",
    SGE_ARCH = "lx-amd64",
    SGE_CELL = "age",
    SGE_ROOT = "/cm/shared/apps/uge/8.8.1"
  )
}

#' Run one workflow with Nextflow and write its run report
#'
#' @return list(status, base, log, queued). Nothing runs when the queue is empty.
#' @noRd
run_pipeline <- function(project, workflow, resume = TRUE, user_asmb = FALSE, echo = TRUE) {
  project <- normalizePath(project)
  workflow <- match.arg(tolower(workflow[1]), c("assemble", "annotate"))
  live <- run_active(project)
  if (!is.na(live)) cli_fail(paste0("A run is still in progress in this project: ", live))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(project, ".sqlite"))
  queued <- tryCatch(length(pipeline_queue(con, workflow, user_asmb)),
                     finally = DBI::dbDisconnect(con))
  if (queued == 0) {
    return(list(status = 0L, base = NA_character_, log = NA_character_, queued = 0L))
  }
  base <- run_basename(workflow)
  log <- file.path(run_dir(project, "nextflow"), paste0(base, ".nextflow.log"))
  args <- nextflow_cmd(workflow, path = project, userAsmbs = user_asmb, base = base)
  if (!resume) args <- stringr::str_remove(args, "-resume")
  res <- processx::run("nextflow", args, wd = project, env = nextflow_env(project),
                       error_on_status = FALSE, echo = echo)
  status <- as.integer(res$status)
  if (file.exists(log)) {
    tryCatch(write_run_report(project, run_report(project, log, exit_status = status, how = "cli")),
             error = function(e) message("Run report not written: ", conditionMessage(e)))
  }
  list(status = status, base = base, log = log, queued = as.integer(queued))
}

#' Write a scheduler job script that runs `command`, and submit it
#'
#' @param command The `mitopilot ...` command line the job runs.
#' @return list(script, lines, submitted, output)
#' @noRd
submit_pipeline <- function(project, command, workflow, dry_run = FALSE) {
  project <- normalizePath(project)
  cfg <- read_config_executor(file.path(project, ".config"))
  if (identical(tolower(cfg$executor), "local")) {
    cli_fail("no scheduler configured; see 'mitopilot cluster setup'")
  }
  base <- run_basename(workflow)
  log <- file.path(run_dir(project, "logs"), paste0(base, ".log"))
  script <- file.path(run_dir(project, "jobs"), paste0(base, ".sh"))
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(script), recursive = TRUE, showWarnings = FALSE)
  lines <- build_submit_script(project, cfg$executor, cfg$queue, command, base, log)
  writeLines(lines, script)
  if (dry_run) return(list(script = script, lines = lines, submitted = FALSE, output = character(0)))
  res <- run_submit(script, cfg$executor, project)
  list(script = script, lines = lines, submitted = isTRUE(res$success), output = res$output)
}

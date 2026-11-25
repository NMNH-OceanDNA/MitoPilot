#' @noRd
cli_parse_cmd <- function(args, cmd, desc, flags = character()) {
  p <- cli_parse_args(args, flags = c(flags, "help"))
  if (isTRUE(p$opts$help)) {
    message(paste0(
      "Usage: mitopilot ", cmd, " [options]\n\nOptions:\n",
      paste0("  --", names(desc), "  ", desc, collapse = "\n")
    ))
    return(NULL)
  }
  cli_check_opts(p$opts, c(names(desc), "help"))
  p
}

#' @noRd
cli_num <- function(x, name) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) cli_fail(paste0("--", name, " must be a number, got '", x, "'"))
  v
}

#' @noRd
cli_cmd_init <- function(args) {
  desc <- c(
    test = "make a test project: reads, assemblies, or annotations",
    mapping = "mapping file (CSV)",
    data = "directory of raw reads",
    assemblies = "directory of user assemblies",
    annotations = "directory of user annotations",
    n = "number of test samples",
    "full-size" = "use full-size test data",
    "no-raw-data" = "project has no raw reads",
    executor = "Nextflow executor (local, slurm, ...)",
    container = "container image",
    force = "overwrite an existing project",
    "no-metadata" = "skip online metadata fetch",
    "min-depth" = "minimum read depth",
    "genetic-code" = "NCBI genetic code table",
    settings = "reserved"
  )
  p <- cli_parse_cmd(args, "init", desc,
    flags = c("full-size", "no-raw-data", "force", "no-metadata"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (!is.null(o[["settings"]])) cli_fail("--settings arrives with 'mitopilot opts'")
  if (length(p$pos) > 1L) cli_fail("init takes at most one directory")

  if (!is.null(o[["test"]])) {
    fn <- switch(o[["test"]],
      reads = "new_test_project",
      assemblies = "new_test_project_userAsmb",
      annotations = "new_test_project_userAnnot",
      cli_fail("--test must be one of: reads, assemblies, annotations")
    )
  } else if (!is.null(o[["annotations"]])) {
    fn <- "new_project_userAnnot"
  } else if (!is.null(o[["assemblies"]])) {
    fn <- "new_project_userAsmb"
  } else {
    fn <- "new_project"
    if (is.null(o[["mapping"]]) || is.null(o[["data"]])) cli_fail("init needs --mapping and --data (or --test)")
  }
  fun <- get(fn, envir = asNamespace("MitoPilot"), mode = "function")
  ok <- setdiff(names(formals(fun)), "...")

  a <- list()
  if (length(p$pos)) a$path <- p$pos[1]
  a$mapping_fn <- o[["mapping"]]
  a$data_path <- o[["data"]]
  a$assembly_path <- o[["assemblies"]]
  a$annotation_path <- o[["annotations"]]
  a$n <- if (!is.null(o[["n"]])) cli_num(o[["n"]], "n")
  a$full_size <- if (isTRUE(o[["full_size"]])) TRUE
  a$no_raw_data <- if (isTRUE(o[["no_raw_data"]])) TRUE
  a$executor <- o[["executor"]]
  a$container <- o[["container"]]
  a$force <- if (isTRUE(o[["force"]])) TRUE
  a$min_depth <- if (!is.null(o[["min_depth"]])) cli_num(o[["min_depth"]], "min-depth")
  a$genetic_code <- if (!is.null(o[["genetic_code"]])) cli_num(o[["genetic_code"]], "genetic-code")
  if (isTRUE(o[["no_metadata"]])) {
    if ("fetch_metadata" %in% ok) {
      a$fetch_metadata <- FALSE
    } else {
      a[c("fetch_geome", "fetch_gbif", "fetch_ncbi")] <- list(FALSE, FALSE, FALSE)
    }
  }
  a$Rproj <- FALSE
  bad <- setdiff(names(a), ok)
  if (length(bad)) {
    cli_fail(paste0("'", fn, "' does not accept: ", paste0("--", gsub("_", "-", bad), collapse = ", ")))
  }
  do.call(fun, a)
  0L
}

#' @noRd
cli_cmd_status <- function(args) {
  desc <- c(
    project = "project directory (default: current directory)",
    json = "print JSON instead of a table",
    stage = "show only this stage (assemble, annotate, export)",
    parked = "show only units waiting for a decision"
  )
  p <- cli_parse_cmd(args, "status", desc, flags = c("json", "parked"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos)) cli_fail(paste0("status takes no arguments, got: ", paste(p$pos, collapse = " ")))
  if (!is.null(o[["stage"]]) && !o[["stage"]] %in% c("assemble", "annotate", "export")) {
    cli_fail("--stage must be one of: assemble, annotate, export")
  }
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  st <- project_status(proj$con)
  parked <- !st$locked & nzchar(st$reasons)
  keep <- rep(TRUE, nrow(st))
  if (!is.null(o[["stage"]])) keep <- keep & st$stage == o[["stage"]]
  if (isTRUE(o[["parked"]])) keep <- keep & parked
  show <- st[keep, , drop = FALSE]
  tally <- warning_tally(proj$con)
  if (isTRUE(o[["json"]])) {
    cat(jsonlite::toJSON(show, dataframe = "rows", na = "null"), "\n", sep = "")
  } else {
    if (nrow(show)) print(show, row.names = FALSE) else cat("No matching units\n")
    message(sprintf("%d units: %d locked, %d parked", nrow(st), sum(st$locked), sum(parked)))
    if (nrow(tally)) {
      message("Warnings: ", paste0(tally$warning, " (", tally$n, ")", collapse = ", "))
    }
  }
  if (any(parked)) 2L else 0L
}

CLI_DESC_PROJECT <- c(project = "project directory (default: current directory)")
CLI_DESC_STAGE <- c(stage = "assemble or annotate (default: inferred per sample)")

#' Units named on the command line
#'
#' A name is a sample ID or `ID:path:scaffold`. Annotate units skip ignored
#' assemblies. `stage` NULL infers it per sample.
#' @return data.frame(ID, path, scaffold, stage, named) with `path`/`scaffold` NA
#'   for assemble units; `named` is TRUE when given in unit spelling.
#' @noRd
cli_check_stage <- function(stage) {
  if (!is.null(stage) && !stage %in% c("assemble", "annotate")) {
    cli_fail("--stage must be one of: assemble, annotate")
  }
}

cli_resolve_units <- function(con, names, stage = NULL) {
  cli_check_stage(stage)
  asm <- DBI::dbGetQuery(con, "SELECT ID, assemble_lock FROM assemble")
  an <- DBI::dbGetQuery(con, paste(
    "SELECT an.ID, an.path, an.scaffold FROM annotate an",
    "LEFT JOIN assemblies s ON s.ID = an.ID AND s.path = an.path AND s.scaffold = an.scaffold",
    "WHERE COALESCE(s.ignore, 0) = 0 ORDER BY an.ID, an.path, an.scaffold"
  ))
  out <- lapply(unique(names), function(nm) {
    parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
    id <- parts[1]
    if (!id %in% asm$ID) cli_fail(paste0("No such sample: ", id))
    u <- an[an$ID == id, , drop = FALSE]
    if (length(parts) == 3L) {
      if (identical(stage, "assemble")) cli_fail(paste0("'", nm, "' names an annotate unit; drop --stage assemble"))
      pa <- suppressWarnings(as.integer(parts[2:3]))
      u <- u[u$path %in% pa[1] & u$scaffold %in% pa[2], , drop = FALSE]
      if (nrow(u) == 0L) cli_fail(paste0("No such unit: ", nm))
      return(data.frame(u, stage = "annotate", named = TRUE))
    }
    if (length(parts) != 1L) cli_fail(paste0("Bad sample or unit name: ", nm, " (use ID or ID:path:scaffold)"))
    st <- stage %||% if (isTRUE(asm$assemble_lock[asm$ID == id] == 1) && nrow(u) > 0L) "annotate" else "assemble"
    if (st == "assemble") {
      return(data.frame(ID = id, path = NA_integer_, scaffold = NA_integer_, stage = "assemble", named = FALSE))
    }
    if (nrow(u) == 0L) cli_fail(paste0("Sample has no annotate units: ", id))
    data.frame(u, stage = "annotate", named = FALSE)
  })
  do.call(rbind, out)
}

cli_unit_label <- function(u) {
  if (is.null(u$path)) return(u$ID)
  ifelse(is.na(u$path), u$ID, paste(u$ID, u$path, u$scaffold, sep = ":"))
}

cli_unit_why <- function(u, why = u$reasons) {
  if (nrow(u) == 0L) return(character())
  paste0(cli_unit_label(u), " (", why, ")")
}

cli_parked_units <- function(con, stage = NULL) {
  cli_check_stage(stage)
  st <- project_status(con)
  st <- st[!st$locked & nzchar(st$reasons), , drop = FALSE]
  if (!is.null(stage)) st <- st[st$stage == stage, , drop = FALSE]
  unique(data.frame(ID = st$ID, path = ifelse(st$stage == "annotate", st$path, NA_integer_),
                    scaffold = ifelse(st$stage == "annotate", st$scaffold, NA_integer_),
                    stage = st$stage, named = FALSE))
}

cli_unit_state_switch <- function(con, u) {
  sw <- if (u$stage[1] == "assemble") {
    DBI::dbGetQuery(con, "SELECT ID, assemble_switch AS sw, assemble_lock AS lk FROM assemble")
  } else {
    DBI::dbGetQuery(con, "SELECT ID, path, scaffold, annotate_switch AS sw, annotate_lock AS lk FROM annotate")
  }
  key <- if (u$stage[1] == "assemble") "ID" else UNIT_KEY
  sw[match(do.call(paste, u[key]), do.call(paste, sw[key])), c("sw", "lk")]
}

cli_by_stage <- function(u, fun) {
  for (st in unique(u$stage)) fun(st, u[u$stage == st, , drop = FALSE])
  invisible(NULL)
}

#' @noRd
cli_cmd_accept <- function(args) {
  desc <- c(CLI_DESC_PROJECT, CLI_DESC_STAGE,
    override = "accept units that fail the pass rules (logged)",
    reason = "note recorded with an override",
    "all-parked" = "every parked unit (needs --override)",
    yes = "do not ask for confirmation")
  p <- cli_parse_cmd(args, "accept SAMPLE...", desc, flags = c("override", "all-parked", "yes"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (!is.null(o[["reason"]]) && !isTRUE(o[["override"]])) cli_fail("--reason needs --override")
  if (isTRUE(o[["all_parked"]])) {
    if (!isTRUE(o[["override"]])) cli_fail("--all-parked needs --override")
    if (length(p$pos)) cli_fail("--all-parked takes no sample names")
  } else if (length(p$pos) == 0L) {
    cli_fail("accept needs sample names or --all-parked")
  }
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  con <- proj$con
  units <- if (isTRUE(o[["all_parked"]])) {
    u <- cli_parked_units(con, o[["stage"]])
    if (nrow(u) == 0L) { message("Nothing parked"); return(0L) }
    u
  } else {
    cli_resolve_units(con, p$pos, o[["stage"]])
  }
  multi <- units[!units$named & units$stage == "annotate", , drop = FALSE]
  multi <- multi[duplicated(multi$ID) | duplicated(multi$ID, fromLast = TRUE), , drop = FALSE]
  if (nrow(multi) && !isTRUE(o[["all_parked"]])) {
    cli_fail(paste0("More than one annotate unit; name one as ID:path:scaffold: ",
                    paste(cli_unit_label(multi), collapse = ", ")))
  }
  settings <- project_settings(con)
  plan <- lapply(split(units, units$stage), function(u) {
    st <- u$stage[1]
    key <- if (st == "assemble") "ID" else UNIT_KEY
    u <- u[key]
    res <- if (st == "assemble") assembly_pass(con, u$ID, settings) else annotation_pass(con, u)
    got <- res[match(do.call(paste, u[key]), do.call(paste, res[key])), , drop = FALSE]
    pass <- got$pass %in% TRUE
    locked <- stage_locked(con, st)
    done <- do.call(paste, u[key]) %in% do.call(paste, locked[key])
    if (any(done)) message("Already accepted: ", paste(cli_unit_label(u[done, , drop = FALSE]), collapse = ", "))
    bad <- u[!pass & !done, , drop = FALSE]
    bad$reasons <- ifelse(is.na(got$reasons), "assembly_not_accepted", got$reasons)[!pass & !done]
    list(stage = st, ok = u[pass & !done, , drop = FALSE], bad = bad)
  })
  n_bad <- sum(vapply(plan, function(x) nrow(x$bad), integer(1)))
  go <- TRUE
  if (n_bad > 0L && isTRUE(o[["override"]])) {
    go <- cli_confirm(sprintf("Override pass rules for %d unit(s)?", n_bad), isTRUE(o[["yes"]]))
  }
  refused <- character()
  declined <- character()
  blocked <- character()
  n_ok <- 0L
  for (x in plan) {
    st <- x$stage
    if (nrow(x$ok)) {
      r <- lock_units(con, proj$dir_out, st, x$ok)
      if (isTRUE(settings$accept_linear)) mark_partial_linear(con, unique(r$written$ID))
      n_ok <- n_ok + nrow(r$written)
      blocked <- c(blocked, cli_unit_why(r$blocked))
    }
    if (nrow(x$bad) == 0L) next
    if (!isTRUE(o[["override"]])) {
      refused <- c(refused, cli_unit_why(x$bad))
    } else if (!go) {
      declined <- c(declined, cli_unit_why(x$bad))
    } else {
      r <- override_accept(con, proj$dir_out, st, x$bad, reason = o[["reason"]] %||% NA)
      n_ok <- n_ok + nrow(r$accepted)
      blocked <- c(blocked, cli_unit_why(r$blocked))
    }
  }
  message(sprintf("Accepted %d unit(s)", n_ok))
  if (length(refused)) message("Refused, fails pass rules (use --override): ", paste(refused, collapse = ", "))
  if (length(declined)) message("Override declined: ", paste(declined, collapse = ", "))
  if (length(blocked)) message("Could not lock: ", paste(blocked, collapse = ", "))
  if (length(refused) || length(declined) || length(blocked)) 2L else 0L
}

#' @noRd
cli_cmd_hold <- function(args) {
  desc <- c(CLI_DESC_PROJECT, CLI_DESC_STAGE)
  p <- cli_parse_cmd(args, "hold SAMPLE...", desc)
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos) == 0L) cli_fail("hold needs sample names")
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  units <- cli_resolve_units(proj$con, p$pos, o[["stage"]])
  cli_by_stage(units, function(st, u) set_unit_state(proj$con, st, u, 0))
  message(sprintf("Held %d unit(s)", nrow(units)))
  0L
}

#' @noRd
cli_cmd_requeue <- function(args) {
  desc <- c(CLI_DESC_PROJECT, CLI_DESC_STAGE, yes = "do not ask for confirmation")
  p <- cli_parse_cmd(args, "requeue SAMPLE...", desc, flags = "yes")
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos) == 0L) cli_fail("requeue needs sample names")
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  con <- proj$con
  units <- cli_resolve_units(con, p$pos, o[["stage"]])
  finished <- 0L
  cli_by_stage(units, function(st, u) {
    s <- cli_unit_state_switch(con, u)
    finished <<- finished + sum(s$sw %in% 2 | s$lk %in% 1)
  })
  if (!cli_confirm(sprintf("%d finished units will rerun", finished), isTRUE(o[["yes"]]))) {
    message("Nothing changed")
    return(1L)
  }
  cli_by_stage(units, function(st, u) {
    locked <- cli_unit_state_switch(con, u)$lk %in% 1
    if (any(locked)) {
      if (st == "assemble") lock_assemblies(con, proj$dir_out, u$ID[locked], lock = FALSE)
      else lock_annotations(con, u[locked, UNIT_KEY], lock = FALSE)
    }
    set_unit_state(con, st, u, 1)
  })
  message(sprintf("Requeued %d unit(s)", nrow(units)))
  0L
}

#' @noRd
cli_cmd_rules <- function(args) {
  desc <- c(CLI_DESC_PROJECT,
    "min-coverage" = "minimum mean read depth to pass",
    "accept-linear" = "let linear assemblies pass (--no-accept-linear to forbid)")
  p <- cli_parse_cmd(args, "rules", desc, flags = "accept-linear")
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos)) cli_fail("rules takes no arguments")
  cov <- if (!is.null(o[["min_coverage"]])) cli_num(o[["min_coverage"]], "min-coverage")
  if (!is.null(cov) && cov < 0) cli_fail("--min-coverage must not be negative")
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  if (!is.null(cov)) set_project_setting(proj$con, "min_mean_depth", cov)
  if (!is.null(o[["accept_linear"]])) set_project_setting(proj$con, "accept_linear", isTRUE(o[["accept_linear"]]))
  s <- project_settings(proj$con)
  message("Pass rules:\n  min-coverage: ", s$min_mean_depth,
          "\n  accept-linear: ", if (s$accept_linear) "yes" else "no")
  0L
}

#' @noRd
cli_cmd_group <- function(args) {
  desc <- c(CLI_DESC_PROJECT,
    samples = "comma-separated samples or ID:path:scaffold units",
    "all-accepted" = "every accepted unit without a group",
    split = "split mixed complete and partial genomes into NAME-complete and NAME-partial")
  p <- cli_parse_cmd(args, "group NAME", desc, flags = c("all-accepted", "split"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos) != 1L) cli_fail("group needs exactly one NAME")
  name <- p$pos[1]
  if (!valid_group_name(name)) {
    cli_fail("Group name may use only letters, numbers, dots, hyphens, and underscores")
  }
  if (!is.null(o[["samples"]]) == isTRUE(o[["all_accepted"]])) {
    cli_fail("group needs exactly one of --samples or --all-accepted")
  }
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  con <- proj$con
  if (isTRUE(o[["all_accepted"]])) {
    u <- DBI::dbGetQuery(con, paste(
      "SELECT an.ID, an.path, an.scaffold FROM annotate an",
      "LEFT JOIN export e ON e.ID = an.ID AND e.path = an.path AND e.scaffold = an.scaffold",
      "WHERE an.annotate_lock = 1 AND e.export_group IS NULL ORDER BY an.ID, an.path, an.scaffold"
    ))
    if (nrow(u) == 0L) { message("No accepted units without a group"); return(0L) }
  } else {
    u <- cli_resolve_units(con, trimws(strsplit(o[["samples"]], ",")[[1]]), "annotate")[UNIT_KEY]
    lk <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate WHERE annotate_lock = 1")
    off <- u[!do.call(paste, u) %in% do.call(paste, lk), , drop = FALSE]
    if (nrow(off)) cli_fail(paste0("Not accepted yet: ", paste(cli_unit_label(off), collapse = ", ")))
  }
  ex <- fetch_export_data(con, session = NULL)
  comp <- ex$completeness[match(do.call(paste, u), do.call(paste, ex[UNIT_KEY]))]
  n_c <- sum(comp == "complete genome", na.rm = TRUE)
  n_p <- sum(comp == "partial genome", na.rm = TRUE)
  groups <- name
  if (n_c > 0L && n_p > 0L) {
    if (!isTRUE(o[["split"]])) {
      cli_fail(sprintf("%d complete and %d partial genomes cannot share a group; use --split", n_c, n_p))
    }
    groups <- ifelse(comp == "complete genome", paste0(name, "-complete"), paste0(name, "-partial"))
  }
  assign_export_group(con, u, groups)
  message(sprintf("%d unit(s) assigned to: %s", nrow(u), paste(unique(groups), collapse = ", ")))
  0L
}

#' @noRd
cli_cmd_metadata <- function(args) {
  desc <- c(CLI_DESC_PROJECT, source = "geome, gbif, or ncbi (default: all)")
  p <- cli_parse_cmd(args, "metadata fetch", desc)
  if (is.null(p)) return(0L)
  o <- p$opts
  if (!identical(p$pos, "fetch")) cli_fail("Usage: mitopilot metadata fetch [--source geome|gbif|ncbi]")
  src <- o[["source"]] %||% c("geome", "gbif", "ncbi")
  if (!all(src %in% c("geome", "gbif", "ncbi"))) cli_fail("--source must be one of: geome, gbif, ncbi")
  proj <- open_project(o[["project"]] %||% ".")
  DBI::dbDisconnect(proj$con)
  for (s in src) {
    message("Fetching ", toupper(s), " metadata")
    fun <- get(paste0("fetch_", s), envir = asNamespace("MitoPilot"), mode = "function")
    fun(path = proj$path)
  }
  0L
}

CLI_JOB_VARS <- c("SLURM_JOB_ID", "JOB_ID", "PBS_JOBID", "LSB_JOBID")

#' The command a scheduler job runs to call mitopilot from a clean shell
#' @noRd
cli_self_command <- function() {
  w <- unname(Sys.which("mitopilot"))
  if (nzchar(w)) return(shQuote(normalizePath(w)))
  paste(shQuote(file.path(R.home("bin"), "Rscript")), "-e",
        shQuote("quit(status = MitoPilot:::cli_main())"))
}

#' @noRd
cli_cmd_run <- function(args) {
  desc <- c(CLI_DESC_PROJECT,
    resume = "resume the previous run (default; --no-resume starts fresh)",
    submit = "submit as a scheduler job instead of running here",
    "dry-run" = "print what would run and stop",
    "auto-accept" = "accept passing units after the run",
    yes = "do not ask for confirmation")
  p <- cli_parse_cmd(args, "run assemble|annotate", desc,
    flags = c("resume", "submit", "dry-run", "auto-accept", "yes"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos) != 1L || !p$pos[1] %in% c("assemble", "annotate")) {
    cli_fail("run needs a workflow: assemble or annotate")
  }
  wf <- p$pos[1]
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  con <- proj$con
  resume <- !identical(o[["resume"]], FALSE)
  dry <- isTRUE(o[["dry_run"]])
  if (!dry) {
    live <- run_active(proj$path)
    if (!is.na(live)) cli_fail(paste0("A run is still in progress in this project: ", live))
  }
  if (wf == "annotate" && annot_import_pending(con)) {
    cli_fail("Imported gene names are not confirmed yet; confirm them in the app first")
  }
  queued <- length(pipeline_queue(con, wf, proj$user_asmb))
  executor <- read_config_executor(file.path(proj$path, ".config"))$executor

  if (isTRUE(o[["submit"]])) {
    cmd <- paste(c(cli_self_command(), "run", wf, if (!resume) "--no-resume",
                   if (isTRUE(o[["auto_accept"]])) "--auto-accept", "--yes",
                   "--project", shQuote(proj$path)), collapse = " ")
    res <- submit_pipeline(proj$path, cmd, wf, dry_run = dry)
    if (dry) {
      cat(res$lines, sep = "\n")
      message(sprintf("Dry run: %d %s unit(s) queued; script %s", queued, wf, res$script))
      return(0L)
    }
    message(if (res$submitted) "Submitted: " else "Submit failed: ", paste(res$output, collapse = " "),
            "\nScript: ", res$script)
    return(if (res$submitted) 0L else 1L)
  }

  if (dry) {
    args <- nextflow_cmd(wf, path = proj$path, userAsmbs = proj$user_asmb)
    if (!resume) args <- stringr::str_remove(args, "-resume")
    cat("nextflow", args, "\n")
    message(sprintf("Dry run: %d %s unit(s) queued", queued, wf))
    return(0L)
  }
  if (!identical(tolower(executor), "local") && !any(nzchar(Sys.getenv(CLI_JOB_VARS)))) {
    message("Warning: executor is '", executor, "' but this is not a scheduler job; ",
            "the Nextflow head process will run on this node (see --submit)")
  }
  status <- 0L
  if (queued == 0L) {
    message("Nothing queued for ", wf)
  } else {
    message(sprintf("Running %s: %d unit(s) queued", wf, queued))
    res <- run_pipeline(proj$path, wf, resume = resume, user_asmb = proj$user_asmb)
    status <- as.integer(res$status)
    where <- if (is.na(res$base)) "" else paste0("; run report: ", run_report_paths(proj$path, res$base)[["txt"]])
    message(if (status == 0L) "Run finished" else paste0("Run failed (exit status ", status, "); log: ", res$log), where)
  }
  if (status != 0L) return(1L)
  if (isTRUE(o[["auto_accept"]])) {
    r <- auto_accept(con, proj$dir_out, wf)
    message(sprintf("Auto-accepted %d unit(s), parked %d", nrow(r$accepted), nrow(r$parked)))
  }
  st <- project_status(con)
  if (any(st$stage == wf & !st$locked & nzchar(st$reasons))) 2L else 0L
}

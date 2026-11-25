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
  if (!is.null(a$executor) && cluster_copy_submit_template(a$executor, a$path %||% ".")) {
    message("Copied the saved '", a$executor, "' submit template into the project")
  }
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
    u <- cli_ungrouped_units(con)
    if (nrow(u) == 0L) { message("No accepted units without a group"); return(0L) }
  } else {
    u <- cli_resolve_units(con, trimws(strsplit(o[["samples"]], ",")[[1]]), "annotate")[UNIT_KEY]
    lk <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate WHERE annotate_lock = 1")
    off <- u[!do.call(paste, u) %in% do.call(paste, lk), , drop = FALSE]
    if (nrow(off)) cli_fail(paste0("Not accepted yet: ", paste(cli_unit_label(off), collapse = ", ")))
  }
  groups <- cli_assign_group(con, u, name, isTRUE(o[["split"]]))
  message(sprintf("%d unit(s) assigned to: %s", nrow(u), paste(unique(groups), collapse = ", ")))
  0L
}

#' Accepted annotation units that have no export group
#' @noRd
cli_ungrouped_units <- function(con) {
  DBI::dbGetQuery(con, paste(
    "SELECT an.ID, an.path, an.scaffold FROM annotate an",
    "LEFT JOIN export e ON e.ID = an.ID AND e.path = an.path AND e.scaffold = an.scaffold",
    "WHERE an.annotate_lock = 1 AND e.export_group IS NULL ORDER BY an.ID, an.path, an.scaffold"
  ))
}

#' Put units in group `name`, or NAME-complete and NAME-partial when genomes mix
#'
#' Members already in the groups are taken into account. Units whose
#' completeness is unknown are left ungrouped (`strict` fails instead), and so
#' are all units when `split` is FALSE and genomes would mix.
#' @return group name per unit, NA for a unit left ungrouped
#' @noRd
cli_assign_group <- function(con, u, name, split, strict = TRUE) {
  ex <- fetch_export_data(con, session = NULL)
  kind <- function(x) c("complete genome" = "complete", "partial genome" = "partial")[x]
  comp <- unname(kind(ex$completeness[match(do.call(paste, u), do.call(paste, ex[UNIT_KEY]))]))
  mem <- function(g) {
    k <- unname(kind(ex$completeness[!is.na(ex$export_group) & ex$export_group == g]))
    k[!is.na(k)]
  }
  groups <- rep(NA_character_, nrow(u))
  unknown <- is.na(comp)
  if (any(unknown) && strict) {
    cli_fail(paste0("Completeness is unknown for: ", paste(cli_unit_label(u[unknown, , drop = FALSE]), collapse = ", ")))
  }
  have <- mem(name)
  batch <- unique(comp[!unknown])
  if (!split) {
    mixed <- unique(c(batch, have))
    if (length(mixed) > 1L) {
      cli_fail(sprintf("%s would mix complete and partial genomes (group holds %d complete, %d partial; adding %d complete, %d partial); use --split",
                       name, sum(have == "complete"), sum(have == "partial"),
                       sum(comp == "complete", na.rm = TRUE), sum(comp == "partial", na.rm = TRUE)))
    }
    groups[!unknown] <- name
  } else {
    sep <- length(c(mem(paste0(name, "-complete")), mem(paste0(name, "-partial")))) > 0L
    for (k in batch) {
      target <- if (length(have)) {
        if (all(have == k) && !sep) name else paste0(name, "-", k)
      } else if (sep || length(batch) > 1L) paste0(name, "-", k) else name
      held <- mem(target)
      if (length(held) && any(held != k)) {
        if (strict) cli_fail(sprintf("Group '%s' already holds %s genomes; cannot add %s genomes", target,
                                     setdiff(held, k)[1], k))
        target <- NA_character_
      }
      groups[!unknown & comp == k] <- target
    }
  }
  ok <- !is.na(groups)
  if (any(ok)) assign_export_group(con, u[ok, , drop = FALSE], groups[ok])
  groups
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

#' Run a workflow when units are queued; returns the run's exit status
#' @noRd
cli_run_queued <- function(proj, wf, queued, resume) {
  if (queued == 0L) {
    message("Nothing queued for ", wf)
    return(0L)
  }
  message(sprintf("Running %s: %d unit(s) queued", wf, queued))
  res <- run_pipeline(proj$path, wf, resume = resume, user_asmb = proj$user_asmb)
  status <- as.integer(res$status)
  base <- res$base %||% NA_character_
  where <- if (is.na(base)) "" else paste0("; run report: ", run_report_paths(proj$path, base)[["txt"]])
  message(if (status == 0L) "Run finished" else paste0("Run failed (exit status ", status, "); log: ", res$log), where)
  status
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
  status <- cli_run_queued(proj, wf, queued, resume)
  if (status != 0L) return(1L)
  if (isTRUE(o[["auto_accept"]])) {
    r <- auto_accept(con, proj$dir_out, wf)
    message(sprintf("Auto-accepted %d unit(s), parked %d", nrow(r$accepted), nrow(r$parked)))
  }
  st <- project_status(con)
  if (any(st$stage == wf & !st$locked & nzchar(st$reasons))) 2L else 0L
}

#' Protein-gene records the outlier check skips: no translation, or the only
#' translated copy of its gene in the group
#' @param a get_export_PCG_annotations() output
#' @noRd
export_unchecked <- function(a) {
  ok <- !is.na(a$translation) & nzchar(a$translation)
  n <- table(a$gene[ok])
  sum(!ok | !a$gene %in% names(n)[n >= 2])
}

CLI_FLAGGED_MENU <- c(
  drop = "drop flagged samples, export the rest",
  keep = "keep all (logged as an override)",
  cancel = "cancel, export nothing"
)

#' @noRd
cli_cmd_export <- function(args) {
  desc <- c(CLI_DESC_PROJECT,
    flagged = "flagged samples: ask, stop, drop, or keep (default: ask in a terminal, stop in a script)",
    yes = "do not ask for confirmation")
  p <- cli_parse_cmd(args, "export GROUP", desc, flags = "yes")
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos) != 1L) cli_fail("export needs exactly one GROUP")
  choice <- o[["flagged"]] %||% if (cli_is_interactive()) "ask" else "stop"
  if (!choice %in% c("ask", "stop", "drop", "keep")) cli_fail("--flagged must be one of: ask, stop, drop, keep")
  if (choice == "ask" && !cli_is_interactive()) cli_fail("--flagged ask needs a terminal; use stop, drop, or keep")
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  cli_export(proj, p$pos[1], choice, isTRUE(o[["yes"]]))
}

#' Export one group of an open project
#' @param choice ask, stop, drop, or keep
#' @param keep_confirmed TRUE when `keep` was already confirmed by the user
#' @noRd
cli_export <- function(proj, group, choice, yes = FALSE, keep_confirmed = FALSE) {
  con <- proj$con
  live <- run_active(proj$path)
  if (!is.na(live)) cli_fail(paste0("A run is still in progress in this project: ", live))

  groups <- DBI::dbGetQuery(con, "SELECT DISTINCT export_group FROM export WHERE export_group IS NOT NULL ORDER BY 1")$export_group
  if (!group %in% groups) {
    cli_fail(paste0("No such export group: ", group,
                    if (length(groups)) paste0(" (groups: ", paste(groups, collapse = ", "), ")")))
  }
  units <- export_group_units(con, group)
  if (nrow(units) == 0L) {
    out <- export_left_out(con, group)
    cli_fail(paste0("Export group '", group, "' has nothing to export; left out: ",
                    paste(sprintf("%s p%s.s%s (%s)", out$ID, out$path, out$scaffold, out$why), collapse = ", ")))
  }
  ids <- unique(units$ID)
  opts <- get_export_opts(con)
  tmpl <- export_last_template(con)
  if (!tmpl %in% list_export_templates(con)) tmpl <- "default"
  if (tmpl != "default") {
    t_opts <- get_export_opts(con, tmpl)
    opts[c("fasta_header", "fasta_header_gene")] <- t_opts[c("fasta_header", "fasta_header_gene")]
  }
  message("Header template: ", tmpl)
  opts$fasta_header <- nmnh_template_pref(con, opts$fasta_header)
  d <- export_table_data(con)
  nmnh <- if (nmnh_template_on(opts$fasta_header)) nmnh_resolve(con, ids)
  bad <- export_header_problem(opts$fasta_header, opts$fasta_header_gene, opts$gene_export,
                               export_with_nmnh(d, nmnh))
  if (!is.null(bad)) cli_fail(paste0("Invalid FASTA header template. ", bad))
  mp <- export_multi_path_ids(d, group)
  if (length(mp)) cli_fail(paste0("Cannot export samples with multiple assembly paths. ", export_multi_path_text(mp)))

  warn <- character()
  if (!is.null(nmnh) && any(!nmnh$ok)) warn <- c(warn, paste("NMNH values missing or wrong.", export_nmnh_text(nmnh)))
  sp <- tryCatch(
    export_specimen_rows(con, ids, c(opts$fasta_header, if (opts$gene_export) opts$fasta_header_gene)),
    error = function(e) { message("Specimen conflict check failed: ", conditionMessage(e)); NULL })
  if (!is.null(sp) && nrow(sp) > 0) {
    warn <- c(warn, paste0("Sample metadata disagrees on items the header template uses:\n",
                           paste(utils::capture.output(print(sp, row.names = FALSE)), collapse = "\n")))
  }
  frag <- export_fragmented_ids(d, group)
  if (length(frag)) warn <- c(warn, paste("Some samples export as multiple records.", export_fragmented_text(frag)))
  if (length(warn)) {
    message(paste0("Warning: ", warn, collapse = "\n\n"))
    if (!cli_confirm("Export anyway?", yes, detail = sub("[.:\n].*", "", warn))) { message("Nothing exported"); return(2L) }
  }

  report <- NULL
  f <- data.frame(ID = character(), gene = character(), why = character())
  if (nrow(units) > 1L) {
    fl <- export_flags(proj, group)
    report <- fl$report
    f <- fl$flagged
    message("Outlier report: ", report)
    n_un <- export_unchecked(get_export_PCG_annotations(con, group))
    if (n_un > 0L) {
      message(sprintf("Note: %d sample-gene pair(s) could not be checked (no translation, or the only copy of the gene in the group)", n_un))
    }
  } else {
    message("Note: one record in the group; the outlier check needs two or more")
  }

  status <- 0L
  export_ids <- NULL
  if (nrow(f) > 0L) {
    flagged_ids <- unique(f$ID)
    message("Flagged by the outlier check:\n", paste(utils::capture.output(print(f, row.names = FALSE)), collapse = "\n"))
    if (choice == "ask") {
      choice <- cli_ask(sprintf("%d sample(s) flagged. What now?", length(flagged_ids)), CLI_FLAGGED_MENU)
    }
    if (choice %in% c("stop", "cancel")) {
      message("Nothing exported", if (choice == "stop") "; choose with --flagged drop or --flagged keep")
      return(2L)
    }
    if (choice == "drop") {
      export_ids <- setdiff(ids, flagged_ids)
      if (length(export_ids) == 0L) {
        message("Nothing to export: every sample in '", group, "' is flagged")
        return(2L)
      }
      message("Dropping flagged: ", paste(flagged_ids, collapse = ", "))
      status <- 2L
    }
    if (choice == "keep") {
      if (!cli_confirm(sprintf("Export %d flagged sample(s) anyway? This is logged as an override.",
                               length(flagged_ids)), yes || keep_confirmed,
                       detail = paste0(length(flagged_ids), " flagged sample(s)"))) {
        message("Nothing exported")
        return(2L)
      }
    }
  }

  kept <- if (!is.null(report) && file.exists(report)) {
    tmp <- tempfile(fileext = ".html")
    file.copy(report, tmp)
    tmp
  }
  export_files(
    group = group, IDs = export_ids,
    fasta_header = opts$fasta_header, fasta_header_gene = opts$fasta_header_gene,
    out_dir = proj$dir_out, generateAAalignments = opts$aa_alignments, gene_export = opts$gene_export,
    review = FALSE, start_aa = opts$start_aa, stop_aa = opts$stop_aa, ident_pct = opts$ident_pct
  )
  if (!is.null(kept)) {
    if (!file.exists(report)) file.copy(kept, report)
    unlink(kept)
  }
  logged <- if (choice == "keep" && nrow(f) > 0L) "override" else if (!is.null(export_ids)) "dropped"
  if (!is.null(logged)) {
    fu <- units[if (logged == "dropped") !units$ID %in% export_ids else units$ID %in% f$ID, , drop = FALSE]
    err <- tryCatch({
      log_review(con, fu, "export", logged,
                 rules = vapply(fu$ID, function(id) pass_collapse(unique(f$why[f$ID == id])), ""))
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(err)) {
      message("Exported, but the ", logged, " was NOT logged for: ",
              paste(cli_unit_label(fu), collapse = ", "), " (", err, ")")
      return(1L)
    }
  }
  if (!is.null(export_ids)) message("Note: the alignment report still shows the dropped samples")
  n <- if (is.null(export_ids)) nrow(units) else sum(units$ID %in% export_ids)
  message(sprintf("Exported %d record(s) to %s", n, file.path(proj$dir_out, "export", group)))
  status
}

#' Groups holding an accepted unit that was never exported and not dropped on purpose
#' @noRd
cli_groups_needing_export <- function(con) {
  g <- DBI::dbGetQuery(con, paste(
    "SELECT DISTINCT e.export_group FROM export e",
    "JOIN annotate an ON an.ID = e.ID AND an.path = e.path AND an.scaffold = e.scaffold",
    "WHERE e.export_group IS NOT NULL AND an.annotate_lock = 1 AND e.export_time_stamp IS NULL",
    "AND COALESCE((SELECT l.action FROM review_log l WHERE l.stage = 'export' AND l.ID = e.ID",
    "AND l.path = e.path AND l.scaffold = e.scaffold ORDER BY l.time_stamp DESC, l.rowid DESC LIMIT 1), '') != 'dropped'",
    "ORDER BY 1"))$export_group
  g[vapply(g, function(x) nrow(export_group_units(con, x)) > 0L, logical(1))]
}

#' @noRd
cli_cmd_auto <- function(args) {
  desc <- c(CLI_DESC_PROJECT,
    group = "export group name (default: the project directory name)",
    flagged = "flagged samples at export: stop, drop, or keep (default: stop)",
    resume = "resume previous runs (default; --no-resume starts fresh)",
    submit = "submit as a scheduler job instead of running here",
    "dry-run" = "print what would run and stop",
    yes = "do not ask for confirmation")
  p <- cli_parse_cmd(args, "auto", desc, flags = c("resume", "submit", "dry-run", "yes"))
  if (is.null(p)) return(0L)
  o <- p$opts
  if (length(p$pos)) cli_fail(paste0("auto takes no arguments, got: ", paste(p$pos, collapse = " ")))
  flagged <- o[["flagged"]] %||% "stop"
  if (!flagged %in% c("stop", "drop", "keep")) cli_fail("--flagged must be one of: stop, drop, keep")
  yes <- isTRUE(o[["yes"]])
  dry <- isTRUE(o[["dry_run"]])
  resume <- !identical(o[["resume"]], FALSE)
  proj <- open_project(o[["project"]] %||% ".")
  on.exit(DBI::dbDisconnect(proj$con), add = TRUE)
  con <- proj$con
  group <- o[["group"]] %||% basename(normalizePath(proj$path))
  if (!valid_group_name(group)) {
    cli_fail(paste0("Group name '", group, "' may use only letters, numbers, dots, hyphens, and underscores; give --group"))
  }
  if (!dry) {
    live <- run_active(proj$path)
    if (!is.na(live)) cli_fail(paste0("A run is still in progress in this project: ", live))
    if (flagged == "keep" && isTRUE(o[["submit"]]) && !yes) {
      cli_fail("--flagged keep with --submit needs --yes")
    }
    if (flagged == "keep") {
      if (!cli_confirm("--flagged keep exports flagged samples and logs an override. Continue?", yes)) {
        message("Cancelled")
        return(1L)
      }
    }
  }
  queued <- c(assemble = length(pipeline_queue(con, "assemble", proj$user_asmb)),
              annotate = length(pipeline_queue(con, "annotate", proj$user_asmb)))

  if (isTRUE(o[["submit"]])) {
    cmd <- paste(c(cli_self_command(), "auto", if (!is.null(o[["group"]])) c("--group", shQuote(group)),
                   if (!is.null(o[["flagged"]])) c("--flagged", flagged), if (!resume) "--no-resume",
                   if (yes) "--yes", "--project", shQuote(proj$path)), collapse = " ")
    res <- submit_pipeline(proj$path, cmd, "auto", dry_run = dry)
    if (dry) {
      cat(res$lines, sep = "\n")
      message(sprintf("Dry run: %d assemble and %d annotate unit(s) queued; script %s",
                      queued[["assemble"]], queued[["annotate"]], res$script))
      return(0L)
    }
    message(if (res$submitted) "Submitted: " else "Submit failed: ", paste(res$output, collapse = " "),
            "\nScript: ", res$script)
    return(if (res$submitted) 0L else 1L)
  }

  if (dry) {
    for (wf in c("assemble", "annotate")) {
      a <- nextflow_cmd(wf, path = proj$path, userAsmbs = proj$user_asmb)
      if (!resume) a <- stringr::str_remove(a, "-resume")
      cat("nextflow", a, "\n")
    }
    message(sprintf("Dry run: %d assemble and %d annotate unit(s) queued; group '%s', flagged: %s",
                    queued[["assemble"]], queued[["annotate"]], group, flagged))
    return(0L)
  }

  executor <- read_config_executor(file.path(proj$path, ".config"))$executor
  if (!identical(tolower(executor), "local") && !any(nzchar(Sys.getenv(CLI_JOB_VARS)))) {
    message("Warning: executor is '", executor, "' but this is not a scheduler job; ",
            "the Nextflow head process will run on this node (see --submit)")
  }
  n_acc <- 0L
  for (wf in c("assemble", "annotate")) {
    if (wf == "annotate" && annot_import_pending(con)) {
      message("Skipping annotate: imported gene names are not confirmed yet; confirm them in the app")
    } else {
      n <- length(pipeline_queue(con, wf, proj$user_asmb))
      if (cli_run_queued(proj, wf, n, resume) != 0L) return(1L)
    }
    r <- auto_accept(con, proj$dir_out, wf)
    n_acc <- n_acc + nrow(r$accepted)
    message(sprintf("Auto-accepted %d %s unit(s), parked %d", nrow(r$accepted), wf, nrow(r$parked)))
  }

  u <- cli_ungrouped_units(con)
  new_groups <- character()
  n_held <- 0L
  if (nrow(u)) {
    g <- cli_assign_group(con, u, group, split = TRUE, strict = FALSE)
    new_groups <- unique(g[!is.na(g)])
    n_held <- sum(is.na(g))
    if (length(new_groups)) {
      message(sprintf("%d unit(s) assigned to: %s", sum(!is.na(g)), paste(new_groups, collapse = ", ")))
    }
    if (n_held) {
      message("Left ungrouped (completeness unknown or would mix genomes): ",
              paste(cli_unit_label(u[is.na(g), , drop = FALSE]), collapse = ", "))
    }
  } else {
    message("No accepted units without a group")
  }
  todo <- cli_groups_needing_export(con)
  results <- integer()
  for (g in todo) {
    message("Exporting group: ", g)
    results[g] <- tryCatch(cli_export(proj, g, flagged, yes, keep_confirmed = flagged == "keep"),
      mitopilot_needs_confirm = function(e) {
        message("Held '", g, "': needs confirmation (", paste(e$detail, collapse = "; "),
                "); run 'mitopilot export ", g, " --yes' or 'mitopilot auto --yes' to proceed")
        2L
      },
      error = function(e) {
        message("Export of '", g, "' failed: ", conditionMessage(e))
        1L
      })
  }

  st <- project_status(con)
  parked <- st[!st$locked & nzchar(st$reasons), , drop = FALSE]
  message(sprintf("Summary: %d unit(s) accepted, %d parked, %d group(s) exported, %d group(s) held back",
                  n_acc, nrow(parked), sum(results == 0L), sum(results != 0L)))
  if (nrow(parked)) {
    codes <- table(unlist(strsplit(parked$reasons, ",")))
    message("Parked: ", paste0(names(codes), " (", as.integer(codes), ")", collapse = ", "))
    message("Run 'mitopilot status --parked' for details.")
  }
  if (any(results == 1L)) 1L else if (nrow(parked) || n_held || any(results == 2L)) 2L else 0L
}

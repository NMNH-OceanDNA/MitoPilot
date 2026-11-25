CLUSTER_STEPS <- c("tools", "config", "submit", "pipeline")

#' @noRd
cluster_ask <- function(prompt, default = NULL) {
  ans <- trimws(cli_readline(paste0(prompt, if (!is.null(default)) paste0(" [", default, "]"), ": ")))
  if (nzchar(ans)) ans else default
}

#' @noRd
cluster_hint <- function(step, name = NULL) {
  nm <- name %||% "NAME"
  switch(step,
    tools = paste0("install or load the missing tools; add 'module load' lines in 'mitopilot cluster edit ", nm, " --submit'"),
    config = paste0("fix the profile with 'mitopilot cluster edit ", nm, "'"),
    submit = paste0("check queue name, account, and time limit in 'mitopilot cluster edit ", nm, " --submit'"),
    pipeline = paste0("check the container engine, cache folder, and queue in 'mitopilot cluster edit ", nm, "'")
  )
}

#' Job script header for a saved profile, tokenized like a project template
#' @noRd
cluster_build_template <- function(scheduler, queue = NULL, account = NULL, extra = NULL,
                                   penv = "mthread", time = "168:00:00", mem = "16G",
                                   setup = character(), native_prefix = NULL) {
  act <- if (!is.null(native_prefix)) file.path(native_prefix, "activate.sh")
  nxf <- if (is.null(act)) nf_pin_version() else tryCatch(native_nf_pin(native_env(act)), error = function(e) NA_character_)
  s <- submission_script(scheduler, queue, "<<NEXTFLOW_CMD>>", "<<JOB_NAME>>", "<<LOG_FILE>>",
                         nxf_ver = nxf, env_setup = act)
  gb <- grepl("^[0-9.]+[Gg]$", mem)
  num <- as.numeric(sub("[A-Za-z]+$", "", mem))
  hm <- if (grepl("^[0-9]+:[0-9]+:[0-9]+$", time)) sub(":[0-9]+$", "", time) else time
  add <- NULL
  if (scheduler == "slurm") {
    s <- sub("^#SBATCH --mem=.*$", paste0("#SBATCH --mem=", mem), s)
    s <- sub("^#SBATCH -t .*$", paste0("#SBATCH -t ", time), s)
    add <- c(if (!is.null(account)) paste0("#SBATCH --account=", account),
             if (!is.null(extra)) paste0("#SBATCH ", extra))
  } else if (scheduler == "pbs") {
    pm <- if (grepl("[A-Za-z]$", mem)) paste0(tolower(mem), if (!grepl("[Bb]$", mem)) "b") else mem
    s <- sub("mem=16gb$", paste0("mem=", pm), s)
    s <- sub("^#PBS -l walltime=.*$", paste0("#PBS -l walltime=", time), s)
    add <- c(if (!is.null(account)) paste0("#PBS -A ", account),
             if (!is.null(extra)) paste0("#PBS ", extra))
  } else if (scheduler == "lsf") {
    mb <- if (gb) format(num * 1000, scientific = FALSE) else mem
    s <- sub("^#BSUB -M .*$", paste0("#BSUB -M ", mb), s)
    s <- sub("^#BSUB -W .*$", paste0("#BSUB -W ", hm), s)
    add <- c(if (!is.null(account)) paste0("#BSUB -P ", account),
             if (!is.null(extra)) paste0("#BSUB ", extra))
  } else if (scheduler == "sge") {
    s <- sub("^#\\$ -pe mthread 1$", paste0("#$ -pe ", penv, " 1"), s)
    add <- c(paste0("#$ -l h_rt=", time), paste0("#$ -l h_vmem=", mem),
             if (!is.null(account)) paste0("#$ -P ", account),
             if (!is.null(extra)) paste0("#$ ", extra))
  }
  if (length(add)) {
    last <- grep("^#(SBATCH|PBS|BSUB|\\$) ", s)
    last <- last[length(last)]
    s <- append(s, add, after = last)
  }
  if (length(setup)) {
    k <- grep("^# (EDIT: load|Native MitoPilot)", s)[1]
    if (!is.na(k)) {
      drop <- if (grepl("EDIT", s[k])) k:(k + 3L) else integer()
      s <- append(if (length(drop)) s[-drop] else s, setup, after = k - 1L)
    }
  }
  s
}

#' Copy a saved profile's submit template into a project that has none
#' @noRd
cluster_copy_submit_template <- function(executor, project) {
  src <- file.path(mitopilot_config_dir(), paste0(executor, ".submit.template"))
  dst <- submit_template_path(project)
  if (!file.exists(src) || file.exists(dst)) return(FALSE)
  file.copy(src, dst)
}

#' @noRd
cluster_context <- function(name = NULL, project = NULL) {
  if (!is.null(name)) {
    cfg <- tryCatch(resolve_config(name), error = function(e) cli_fail(conditionMessage(e)))
    lines <- readLines(cfg)
    tmpl <- file.path(dirname(cfg), paste0(name, ".submit.template"))
    base <- mitopilot_config_dir()
    project <- NULL
  } else {
    project <- normalizePath(project %||% ".")
    cfg <- file.path(project, ".config")
    if (!file.exists(cfg)) cli_fail(paste0("No .config in ", project, "; give a profile NAME or --project"))
    lines <- readLines(cfg)
    tmpl <- submit_template_path(project)
    base <- project
  }
  lines <- fill_config(lines, list(
    CONTAINER_ID = paste0("macguigand/mitopilot:", utils::packageVersion("MitoPilot")),
    RAW_DIR = "NA", ASMB_DIR = "NA", MIN_DEPTH = "2000000", NCBI_API_KEY = ""
  ))
  scratch <- file.path(base, paste0(".mitopilot_check_", Sys.getpid()))
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  conf <- file.path(scratch, "nextflow.config")
  writeLines(lines, conf)
  ex <- read_config_executor(conf)
  eng <- sub("^\\s*([a-z]+).*$", "\\1", extract_container_engine(lines)[1])
  list(name = name, project = project, config = conf, scratch = scratch,
       template = if (file.exists(tmpl)) tmpl, executor = ex$executor, queue = ex$queue,
       native = ex$native_activate, engine = eng, lines = lines)
}

#' @noRd
cluster_env <- function(ctx) {
  nat <- if (!is.null(ctx$native) && file.exists(ctx$native)) tryCatch(native_env(ctx$native), error = function(e) NULL)
  pin <- if (length(nat)) native_nf_pin(nat) else nf_pin_version()
  c("current", nat, if (!is.na(pin)) c(NXF_VER = pin), NXF_ANSI_SUMMARY = TRUE, COLUMNS = 500)
}

#' @noRd
cluster_path_from_native <- function(ctx) {
  old <- Sys.getenv("PATH")
  if (!is.null(ctx$native) && file.exists(ctx$native)) {
    nat <- tryCatch(native_env(ctx$native), error = function(e) NULL)
    if (!is.null(nat)) Sys.setenv(PATH = nat[["PATH"]])
  }
  old
}

#' @noRd
cluster_tools_found <- function(ctx) {
  old <- cluster_path_from_native(ctx)
  on.exit(Sys.setenv(PATH = old), add = TRUE)
  tools <- c(submit_command(ctx$executor), "nextflow", "java", if (is.null(ctx$native)) ctx$engine)
  tools <- setdiff(tools, "none")
  found <- vapply(tools, function(x) nzchar(Sys.which(x)), logical(1))
  if (!is.null(ctx$native)) found[["activate.sh"]] <- file.exists(ctx$native)
  found
}

#' @noRd
cluster_nextflow_version <- function(ctx) {
  old <- cluster_path_from_native(ctx)
  on.exit(Sys.setenv(PATH = old), add = TRUE)
  v <- nf_installed_version()
  st <- nf_version_status(v)
  if (st == "missing") return(list(ok = FALSE, detail = "'nextflow -version' reported no version; check Java 17+ is on PATH"))
  if (st == "too_old") return(list(ok = FALSE, detail = paste0("Nextflow ", v, " is too old. ", nf_update_hint())))
  list(ok = TRUE, detail = paste0("nextflow ", v))
}

#' @noRd
cluster_tail <- function(x, n = 15L) {
  x <- strsplit(trimws(x), "\n", fixed = TRUE)[[1]]
  paste(utils::tail(x, n), collapse = "\n")
}

#' @noRd
cluster_config_parses <- function(ctx) {
  res <- processx::run("nextflow", c("config", "."), wd = ctx$scratch, env = cluster_env(ctx),
                       error_on_status = FALSE)
  list(ok = identical(res$status, 0L), detail = if (res$status == 0L) "" else cluster_tail(paste(res$stderr, res$stdout)))
}

#' @noRd
cluster_submit_lines <- function(ctx, command, job, log) {
  if (!is.null(ctx$template)) return(fill_submit_template(readLines(ctx$template), command, job, log))
  if (!is.null(ctx$project)) return(build_submit_script(ctx$project, ctx$executor, ctx$queue, command, job, log))
  submission_script(ctx$executor, ctx$queue, command, job, log, env_setup = ctx$native)
}

#' @noRd
cluster_submit_job <- function(ctx, wait) {
  if (is.null(submit_command(ctx$executor))) {
    return(list(ok = TRUE, detail = paste0("executor '", ctx$executor, "' has no scheduler; skipped")))
  }
  job <- "mitopilot_check"
  log <- file.path(ctx$scratch, "check.log")
  script <- file.path(ctx$scratch, "check.sh")
  writeLines(cluster_submit_lines(ctx, "hostname", job, log), script)
  res <- run_submit(script, ctx$executor, ctx$scratch)
  if (!isTRUE(res$success)) return(list(ok = FALSE, detail = res$output))
  deadline <- Sys.time() + wait
  repeat {
    txt <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
    hit <- grep("nextflow exited with status", txt, value = TRUE)
    if (length(hit)) {
      ok <- grepl("status 0$", hit[1])
      return(list(ok = ok, detail = if (ok) paste0("job ran on ", txt[nzchar(txt) & !grepl("^(---|\\[MitoPilot\\])", txt)][1])
                                    else paste(txt, collapse = "\n")))
    }
    if (Sys.time() > deadline) {
      return(list(ok = FALSE, detail = paste0("no result in ", wait, " seconds; scheduler said: ",
                                              paste(res$output, collapse = " "))))
    }
    Sys.sleep(2)
  }
}

#' @noRd
cluster_pipeline_runs <- function(ctx, wait = 1800) {
  wd <- file.path(ctx$scratch, "run")
  dir.create(wd, showWarnings = FALSE)
  res <- processx::run("nextflow",
    c("run", app_sys("cli", "check.nf"), "-c", ctx$config, "-w", file.path(wd, "work"), "-ansi-log", "false"),
    wd = wd, env = cluster_env(ctx), error_on_status = FALSE, timeout = wait)
  list(ok = identical(res$status, 0L),
       detail = if (res$status == 0L) "MitoPilot loads on a compute node" else cluster_tail(paste(res$stderr, res$stdout), 25L))
}

#' @noRd
cluster_check <- function(name = NULL, project = NULL, wait = 300) {
  ctx <- cluster_context(name, project)
  on.exit(unlink(ctx$scratch, recursive = TRUE), add = TRUE)
  steps <- list(
    tools = function() {
      found <- cluster_tools_found(ctx)
      if (!all(found)) return(list(ok = FALSE, detail = paste("not found:", paste(names(found)[!found], collapse = ", "))))
      v <- cluster_nextflow_version(ctx)
      if (v$ok) v$detail <- paste("found:", paste(names(found), collapse = ", "), "(", v$detail, ")")
      v
    },
    config = function() cluster_config_parses(ctx),
    submit = function() cluster_submit_job(ctx, wait),
    pipeline = function() cluster_pipeline_runs(ctx)
  )
  out <- data.frame(step = character(), ok = logical(), detail = character())
  for (s in names(steps)) {
    r <- steps[[s]]()
    out[nrow(out) + 1L, ] <- list(s, isTRUE(r$ok), paste(r$detail, collapse = "\n"))
    if (!isTRUE(r$ok)) break
  }
  out
}

#' @noRd
cluster_saved_path <- function(name) {
  p <- file.path(mitopilot_config_dir(), paste0("config.", name))
  if (!file.exists(p)) cli_fail(paste0("No saved profile '", name, "'; see 'mitopilot cluster list'"))
  p
}

#' @noRd
cluster_setup <- function(name = NULL) {
  if (!cli_is_interactive()) cli_fail("cluster setup needs a terminal")
  dir <- mitopilot_config_dir()
  name <- name %||% cluster_ask("Profile name")
  if (is.null(name) || !grepl("^[A-Za-z0-9_.-]+$", name)) cli_fail("Profile name must be letters, digits, '_', '-', or '.'")
  overwrite <- file.exists(file.path(dir, paste0("config.", name)))
  if (overwrite && !cli_confirm(paste0("Profile '", name, "' exists; replace it?"))) return(1L)
  sched <- cli_ask("Scheduler:", c(slurm = "slurm", sge = "sge", pbs = "pbs", lsf = "lsf"))
  queue <- cluster_ask("Queue or partition (blank to skip)")
  account <- cluster_ask("Account or project code (blank to skip)")
  extra <- cluster_ask("Extra scheduler options (blank to skip)")
  penv <- if (sched == "sge") cluster_ask("Parallel environment", "mthread") else "mthread"
  engine <- cli_ask("Container engine:", c(singularity = "singularity", apptainer = "apptainer",
                                           docker = "docker", none = "none (native install)"))
  folder <- cluster_ask(if (engine == "none") "Native install folder" else "Container cache folder (blank to skip)")
  if (engine == "none" && is.null(folder)) cli_fail("A native install needs its folder")
  time <- cluster_ask("Head-job time limit", "168:00:00")
  mem <- cluster_ask("Head-job memory", "16G")
  message("Setup lines to run before the pipeline, one per line (blank line to finish):")
  setup <- character()
  repeat {
    l <- trimws(cli_readline("> "))
    if (!nzchar(l)) break
    setup <- c(setup, l)
  }
  path <- generate_config(name, scheduler = sched, container_engine = engine,
    container_cache = if (engine != "none") folder, native_prefix = if (engine == "none") folder,
    queue = queue, account = account, cluster_options = extra, penv = penv,
    profile_dir = dir, overwrite = overwrite)
  tmpl <- cluster_build_template(sched, queue, account, extra, penv, time, mem, setup,
                                 native_prefix = if (engine == "none") folder)
  writeLines(tmpl, file.path(dir, paste0(name, ".submit.template")))
  message("Saved submit template to ", file.path(dir, paste0(name, ".submit.template")))
  if (cli_confirm("Run 'cluster check' now?")) return(cluster_report(name, NULL, 300))
  0L
}

#' @noRd
cluster_report <- function(name, project, wait) {
  r <- cluster_check(name, project, wait)
  for (i in seq_len(nrow(r))) {
    message(sprintf("%-9s %s", r$step[i], if (r$ok[i]) "ok" else "FAILED"),
            if (nzchar(r$detail[i])) paste0("\n  ", gsub("\n", "\n  ", r$detail[i])))
  }
  if (all(r$ok) && nrow(r) == length(CLUSTER_STEPS)) {
    message("All checks passed")
    return(0L)
  }
  message("Hint: ", cluster_hint(r$step[nrow(r)], name))
  1L
}

#' @noRd
cli_cmd_cluster <- function(args) {
  desc <- c(
    name = "profile name (setup)",
    project = "project directory (check, apply)",
    submit = "edit the submit template instead of the config",
    yes = "do not ask for confirmation",
    wait = "seconds to wait for the test job (default 300)"
  )
  p <- cli_parse_cmd(args, "cluster setup|list|show NAME|edit NAME|apply NAME|check [NAME]", desc,
                     flags = c("submit", "yes"))
  if (is.null(p)) return(0L)
  o <- p$opts
  sub <- p$pos[1]
  nm <- p$pos[2]
  if (length(p$pos) > 2L) cli_fail("Too many arguments")
  need <- function() if (is.na(nm)) cli_fail(paste0("cluster ", sub, " needs a profile NAME"))
  switch(sub %||% "",
    setup = cluster_setup(o[["name"]]),
    list = {
      d <- list_configs()
      d$submit_template <- file.exists(file.path(dirname(d$path), paste0(d$name, ".submit.template")))
      cat(utils::capture.output(print(d[c("name", "type", "submit_template", "path")], row.names = FALSE)), sep = "\n")
      0L
    },
    show = {
      need()
      cfg <- tryCatch(resolve_config(nm), error = function(e) cli_fail(conditionMessage(e)))
      cat("== ", cfg, " ==\n", sep = "")
      cat(readLines(cfg), sep = "\n")
      t <- file.path(dirname(cfg), paste0(nm, ".submit.template"))
      if (file.exists(t)) {
        cat("\n== ", t, " ==\n", sep = "")
        cat(readLines(t), sep = "\n")
      }
      0L
    },
    edit = {
      need()
      if (!cli_is_interactive()) cli_fail("cluster edit needs a terminal")
      f <- if (isTRUE(o[["submit"]])) file.path(mitopilot_config_dir(), paste0(nm, ".submit.template")) else cluster_saved_path(nm)
      if (!file.exists(f)) cli_fail(paste0("No such file: ", f))
      utils::file.edit(f, editor = Sys.getenv("EDITOR", "vi"))
      0L
    },
    apply = {
      need()
      cluster_saved_path(nm)
      proj <- normalizePath(o[["project"]] %||% ".")
      if (!file.exists(file.path(proj, ".config"))) cli_fail(paste0("No .config in ", proj))
      if (!cli_confirm(paste0("Replace the .config and submit template in ", proj, " with profile '", nm, "'?"), isTRUE(o[["yes"]]))) {
        message("Nothing changed")
        return(1L)
      }
      if (!isTRUE(migrate_config(proj, nm))) cli_fail("Could not regenerate .config")
      tp <- submit_template_path(proj)
      if (file.exists(tp)) unlink(tp)
      if (cluster_copy_submit_template(nm, proj)) message("Copied submit template")
      0L
    },
    check = cluster_report(if (is.na(nm)) NULL else nm, o[["project"]],
                           if (is.null(o[["wait"]])) 300 else cli_num(o[["wait"]], "wait")),
    cli_fail("cluster needs one of: setup, list, show, edit, apply, check")
  )
}

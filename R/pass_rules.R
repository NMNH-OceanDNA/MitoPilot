PASS_ASSEMBLE_POOR_REF <- c("poor", "failed")

#' Stored per-base depth string to a numeric vector
#'
#' `assemblies.depth` holds one value per base, space separated.
#' @return Numeric vector, empty when the value is missing or blank.
#' @noRd
parse_stored_depth <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(numeric())
  suppressWarnings(as.numeric(strsplit(x, " ")[[1]]))
}

pass_collapse <- function(codes) paste(codes, collapse = ",")

#' Does each sample's non-ignored assembly units have enough depth?
#' @noRd
depth_below <- function(depth, min_depth) {
  if (length(depth) == 0L) return(TRUE)
  any(vapply(depth, function(d) {
    v <- parse_stored_depth(d)
    length(v) == 0L || anyNA(v) || !isTRUE(mean(v) >= min_depth)
  }, logical(1)))
}

#' Which assembled samples are clean enough to accept without review
#'
#' @return data.frame(ID, pass, reasons); `reasons` is a comma-separated string
#'   of codes, empty when passing.
#' @noRd
assembly_pass <- function(con, ids = NULL, settings = project_settings(con)) {
  asm <- DBI::dbGetQuery(con, "SELECT ID, assemble_switch, topology, poor_blast_ref FROM assemble ORDER BY ID")
  units <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold, depth FROM assemblies WHERE COALESCE(ignore, 0) = 0")
  if (!is.null(ids)) asm <- asm[asm$ID %in% ids, , drop = FALSE]
  reasons <- vapply(asm$ID, function(id) {
    i <- match(id, asm$ID)
    u <- units[units$ID == id, , drop = FALSE]
    sw <- asm$assemble_switch[i]
    codes <- c(
      not_run = is.na(sw) || !sw %in% c(2, 3),
      failed = !is.na(sw) && sw == 3,
      multiple_paths = length(unique(u$path)) > 1L,
      multiple_scaffolds = length(unique(u$scaffold)) > 1L,
      linear = !isTRUE(settings$accept_linear) && !identical(tolower(asm$topology[i]), "circular"),
      low_coverage = depth_below(u$depth, settings$min_mean_depth),
      poor_reference = isTRUE(asm$poor_blast_ref[i] %in% PASS_ASSEMBLE_POOR_REF)
    )
    pass_collapse(names(codes)[codes])
  }, character(1))
  data.frame(ID = asm$ID, pass = !nzchar(reasons), reasons = unname(reasons))
}

#' Which assembly-locked annotate units are clean enough to accept
#'
#' @return data.frame(ID, path, scaffold, pass, reasons)
#' @noRd
annotation_pass <- function(con, units = NULL) {
  an <- DBI::dbGetQuery(con, paste(
    "SELECT an.ID, an.path, an.scaffold, an.annotate_switch, an.warnings, an.problematic",
    "FROM annotate an JOIN assemble a ON a.ID = an.ID",
    "LEFT JOIN assemblies s ON s.ID = an.ID AND s.path = an.path AND s.scaffold = an.scaffold",
    "WHERE a.assemble_lock = 1 AND COALESCE(s.ignore, 0) = 0",
    "ORDER BY an.ID, an.path, an.scaffold"
  ))
  if (!is.null(units)) {
    an <- an[do.call(paste, an[UNIT_KEY]) %in% do.call(paste, units[UNIT_KEY]), , drop = FALSE]
  }
  reasons <- vapply(seq_len(nrow(an)), function(i) {
    sw <- an$annotate_switch[i]
    w <- an$warnings[i]
    codes <- c(
      not_run = is.na(sw) || !sw %in% c(2, 3),
      failed = !is.na(sw) && sw == 3,
      warnings = is.na(w) || w > 0,
      problematic = identical(an$problematic[i], "yes")
    )
    pass_collapse(names(codes)[codes])
  }, character(1))
  data.frame(an[UNIT_KEY], pass = !nzchar(reasons), reasons = unname(reasons))
}

#' Outlier flags of an export group, plus the HTML alignment report
#'
#' @return list(flagged = data.frame(ID, gene, why), result, report)
#' @noRd
export_flags <- function(project, group) {
  con <- project$con
  o <- get_export_opts(con)
  out_path <- file.path(project$dir_out, "export", group)
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  make_PCG_alignments(group, project$db, out_path, o$start_aa, o$stop_aa, o$ident_pct)
  res <- flag_PCG_outliers(group, project$db, o$start_aa, o$stop_aa, o$ident_pct)
  f <- res$flags
  kinds <- list(start_offset = f$start_flag, stop_offset = f$stop_flag,
                low_identity = f$identity_flag, internal_stop = f$internal_stop_flag)
  flagged <- do.call(rbind, c(
    lapply(names(kinds), function(k) data.frame(ID = f$ID[kinds[[k]]], gene = f$gene[kinds[[k]]],
                                                why = rep(k, sum(kinds[[k]])))),
    list(data.frame(ID = res$internal_stops$ID, gene = res$internal_stops$gene,
                    why = rep("internal_stop", nrow(res$internal_stops))))
  ))
  flagged <- unique(flagged)
  flagged <- flagged[order(flagged$ID, flagged$gene, flagged$why), , drop = FALSE]
  rownames(flagged) <- NULL
  list(flagged = flagged, result = res,
       report = file.path(out_path, paste0("AA_alignments_", group, ".html")))
}

#' Record accepted units in the review log
#' @param units data.frame with `ID`, and `path`, `scaffold` for annotate units.
#' @param rules One string or one per unit.
#' @noRd
log_review <- function(con, units, stage, action = c("auto", "override", "dropped"), rules = "", reason = NA) {
  action <- match.arg(action)
  n <- nrow(units)
  if (n == 0L) return(invisible(NULL))
  na_col <- function(nm) if (nm %in% names(units)) units[[nm]] else rep(NA_integer_, n)
  DBI::dbAppendTable(con, "review_log", data.frame(
    ID = units$ID, path = na_col("path"), scaffold = na_col("scaffold"),
    stage = stage, action = action, rules = rep_len(as.character(rules), n),
    reason = rep_len(as.character(reason), n), time_stamp = as.integer(Sys.time())
  ))
  invisible(NULL)
}

mark_partial_linear <- function(con, ids) {
  for (id in ids) {
    DBI::dbExecute(con, paste(
      "UPDATE annotate SET partial = 'yes' WHERE ID = ? AND",
      "COALESCE((SELECT topology FROM assemble WHERE ID = ?), '') != 'circular'"
    ), params = list(id, id))
  }
}

#' Lock the units of one stage and report what went through
#'
#' Annotate units are locked one sample at a time so one sample's multi-path
#' block does not stop the others.
#' @return list(written = units locked, blocked = data.frame(ID, [path, scaffold], reasons))
#' @noRd
lock_units <- function(con, dir_out, stage, units) {
  if (stage == "assemble") {
    r <- lock_assemblies(con, dir_out, units$ID)
    blocked <- units[units$ID %in% r$stale$ID, , drop = FALSE]
    written <- units[!units$ID %in% r$stale$ID, , drop = FALSE]
    blocked$reasons <- rep("output_missing", nrow(blocked))
  } else {
    written <- units[0, , drop = FALSE]
    blocked <- units[0, , drop = FALSE]
    for (id in unique(units$ID)) {
      u <- units[units$ID == id, , drop = FALSE]
      r <- lock_annotations(con, u)
      if (length(r$multi_path) > 0L) {
        u$reasons <- rep("multiple_paths", nrow(u))
        blocked <- rbind(blocked, u)
      } else {
        written <- rbind(written, u)
      }
    }
    if (nrow(blocked) == 0L) blocked$reasons <- character()
  }
  list(written = written, blocked = blocked)
}

stage_locked <- function(con, stage) {
  if (stage == "assemble") {
    DBI::dbGetQuery(con, "SELECT ID FROM assemble WHERE assemble_lock = 1")
  } else {
    DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate WHERE annotate_lock = 1")
  }
}

#' Lock units and log them; a failed log write undoes the lock
#' @return list(written, blocked)
#' @noRd
lock_and_log <- function(con, dir_out, stage, units, action, rules = NULL, reason = NA, partial = FALSE) {
  key <- if (stage == "assemble") "ID" else UNIT_KEY
  r <- lock_units(con, dir_out, stage, units)
  if (nrow(r$written) == 0L) return(r)
  tryCatch({
    if (partial) mark_partial_linear(con, unique(r$written$ID))
    rl <- if (is.null(rules)) "" else {
      x <- rules$reasons[match(do.call(paste, r$written[key]), do.call(paste, rules[key]))]
      ifelse(is.na(x), "", x)
    }
    log_review(con, r$written, stage, action, rules = rl, reason = reason)
  }, error = function(e) {
    if (stage == "assemble") lock_assemblies(con, dir_out, r$written$ID, lock = FALSE)
    else lock_annotations(con, r$written, lock = FALSE)
    stop(e)
  })
  r
}

#' Lock every passing, not yet locked unit and log it as auto
#'
#' @return list(accepted, parked); parked units carry their reason codes.
#' @noRd
auto_accept <- function(con, dir_out, stage = c("assemble", "annotate")) {
  stage <- match.arg(stage)
  settings <- project_settings(con)
  res <- if (stage == "assemble") assembly_pass(con, settings = settings) else annotation_pass(con)
  key <- if (stage == "assemble") "ID" else UNIT_KEY
  locked <- stage_locked(con, stage)
  res <- res[!do.call(paste, res[key]) %in% do.call(paste, locked[key]), , drop = FALSE]
  cand <- res[res$pass, key, drop = FALSE]
  parked <- res[!res$pass, c(key, "reasons"), drop = FALSE]
  if (stage == "assemble" && !dir.exists(dir_out)) {
    blocked <- cand
    blocked$reasons <- rep("output_missing", nrow(blocked))
    r <- list(written = cand[0, , drop = FALSE], blocked = blocked)
  } else {
    r <- lock_and_log(con, dir_out, stage, cand, "auto", partial = isTRUE(settings$accept_linear))
  }
  parked <- rbind(parked, r$blocked[c(key, "reasons")])
  rownames(parked) <- NULL
  list(accepted = r$written, parked = parked)
}

#' Lock units regardless of the pass rules and log the failing codes
#'
#' Units that do not exist, or whose sample is not assembly-locked, come back
#' in `blocked` and are neither locked nor logged.
#' @return list(accepted, blocked)
#' @noRd
override_accept <- function(con, dir_out, stage = c("assemble", "annotate"), units, reason = NA) {
  stage <- match.arg(stage)
  key <- if (stage == "assemble") "ID" else UNIT_KEY
  units <- unique(units[key])
  res <- if (stage == "assemble") assembly_pass(con) else annotation_pass(con)
  known <- if (stage == "assemble") {
    DBI::dbGetQuery(con, "SELECT ID FROM assemble")
  } else {
    DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate")
  }
  k <- do.call(paste, units[key])
  unknown <- !k %in% do.call(paste, known[key])
  not_ok <- !unknown & !k %in% do.call(paste, res[key])
  pre <- rbind(
    cbind(units[unknown, , drop = FALSE], reasons = rep("unknown_unit", sum(unknown))),
    cbind(units[not_ok, , drop = FALSE], reasons = rep("assembly_not_accepted", sum(not_ok)))
  )
  ok <- units[!unknown & !not_ok, , drop = FALSE]
  r <- lock_and_log(con, dir_out, stage, ok, "override", rules = res, reason = reason)
  blocked <- rbind(pre, r$blocked)
  rownames(blocked) <- NULL
  list(accepted = r$written, blocked = blocked)
}

#' One row per unit with its furthest stage, lock, and who accepted it
#' @noRd
project_status <- function(con) {
  asm <- DBI::dbGetQuery(con, "SELECT ID, assemble_switch, assemble_lock FROM assemble ORDER BY ID")
  an <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold, annotate_switch, annotate_lock FROM annotate ORDER BY ID, path, scaffold")
  ex <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold, export_group FROM export")
  log <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold, stage, action FROM review_log ORDER BY time_stamp, rowid")
  ap <- assembly_pass(con)
  np <- annotation_pass(con)
  label <- function(code) {
    vapply(code, function(k) if (is.na(k)) NA_character_ else MP_STATE_META[[as.character(k)]]$label %||% NA_character_, character(1))
  }
  accepted_by <- function(stage, id, path = NA, scaffold = NA) {
    l <- log[log$stage == stage & log$ID == id, , drop = FALSE]
    if (stage == "annotate") l <- l[l$path %in% path & l$scaffold %in% scaffold, , drop = FALSE]
    if (nrow(l) == 0L) "manual" else l$action[nrow(l)]
  }
  rows <- lapply(seq_len(nrow(asm)), function(i) {
    id <- asm$ID[i]
    a <- an[an$ID == id, , drop = FALSE]
    if (nrow(a) == 0L) a <- data.frame(ID = id, path = NA_integer_, scaffold = NA_integer_,
                                       annotate_switch = NA_integer_, annotate_lock = NA_integer_)
    do.call(rbind, lapply(seq_len(nrow(a)), function(j) {
      p <- a$path[j]; s <- a$scaffold[j]
      g <- ex$export_group[ex$ID == id & ex$path %in% p & ex$scaffold %in% s]
      g <- if (length(g) > 0L) g[1] else NA_character_
      asm_locked <- isTRUE(asm$assemble_lock[i] == 1)
      stage <- if (!is.na(g)) "export" else if (asm_locked) "annotate" else "assemble"
      if (stage == "assemble" || is.na(p)) {
        sw <- asm$assemble_switch[i]; lk <- asm_locked
        why <- if (lk) "" else ap$reasons[ap$ID == id]
        who_stage <- "assemble"
      } else {
        sw <- a$annotate_switch[j]; lk <- isTRUE(a$annotate_lock[j] == 1)
        why <- if (stage == "annotate" && !lk && !is.na(p)) {
          np$reasons[np$ID == id & np$path == p & np$scaffold == s]
        } else ""
        who_stage <- "annotate"
      }
      if (length(why) == 0L) why <- ""
      data.frame(ID = id, path = p, scaffold = s, stage = stage, state = label(sw),
                 locked = lk, accepted_by = if (lk) accepted_by(who_stage, id, p, s) else NA_character_,
                 reasons = why, export_group = g)
    }))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) out <- data.frame(ID = character(), path = integer(), scaffold = integer(),
    stage = character(), state = character(), locked = logical(), accepted_by = character(),
    reasons = character(), export_group = character())
  rownames(out) <- NULL
  out
}

#' Count of each curation warning across all annotations
#' @return data.frame(warning, n)
#' @noRd
warning_tally <- function(con) {
  w <- DBI::dbGetQuery(con, "SELECT warnings FROM annotations WHERE warnings IS NOT NULL AND warnings != ''")$warnings
  w <- trimws(unlist(strsplit(w, ";", fixed = TRUE)))
  w <- w[nzchar(w)]
  if (length(w) == 0L) return(data.frame(warning = character(), n = integer()))
  t <- sort(table(w), decreasing = TRUE)
  data.frame(warning = names(t), n = as.integer(t))
}

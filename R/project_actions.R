UNIT_KEY <- c("ID", "path", "scaffold")

#' Set the state of assemble samples or annotate units
#'
#' @param con An open connection to a project database.
#' @param stage "assemble" or "annotate".
#' @param units Data frame of units; `ID` for assemble, `ID`, `path`,
#'   `scaffold` for annotate.
#' @param state State code, one of 0:4.
#'
#' @return The data frame of rows written.
#' @noRd
set_unit_state <- function(con, stage = c("assemble", "annotate"), units, state) {
  stage <- match.arg(stage)
  stopifnot(length(state) == 1, state %in% 0:4)
  by <- if (stage == "assemble") "ID" else UNIT_KEY
  upd <- units[by]
  upd[[paste0(stage, "_switch")]] <- rep(as.numeric(state), nrow(upd))
  dplyr::tbl(con, stage) |>
    dplyr::rows_update(upd, unmatched = "ignore", in_place = TRUE, copy = TRUE, by = by)
  upd
}

#' Lock or unlock assemble samples
#'
#' Locking holds back samples whose published assembly output is not on disk
#' and clears any queued scaffold join.
#'
#' @return list(written, stale): the rows written, and the
#'   [stale_assemble_dirs()] rows of samples that were not locked.
#' @noRd
lock_assemblies <- function(con, dir_out, ids, lock = TRUE) {
  stale <- if (lock) {
    tryCatch(
      stale_assemble_dirs(con, dir_out, ids = ids, pending_only = FALSE),
      error = function(e) NULL
    )
  }
  if (is.null(stale)) stale <- data.frame(ID = character(0))
  keep <- ids[!ids %in% stale$ID]
  upd <- data.frame(ID = keep, assemble_lock = rep(as.numeric(lock), length(keep)))
  # A locked sample is never admitted by WF1, so a queued join could never run.
  if (lock) upd$join_switch <- rep(NA_integer_, nrow(upd))
  upd$time_stamp <- rep(as.integer(Sys.time()), nrow(upd))
  if (nrow(upd) > 0L) {
    dplyr::tbl(con, "assemble") |>
      dplyr::rows_update(upd, unmatched = "ignore", in_place = TRUE, copy = TRUE, by = "ID")
  }
  list(written = upd, stale = stale)
}

#' Lock or unlock annotate units
#'
#' A sample can export only one assembly path, so a lock that would leave any
#' sample with more than one path locked writes nothing.
#'
#' @return list(written, multi_path): the rows written, and the IDs that
#'   blocked the lock.
#' @noRd
lock_annotations <- function(con, units, lock = TRUE) {
  upd <- units[UNIT_KEY]
  upd$annotate_lock <- rep(as.numeric(lock), nrow(upd))
  if (lock) {
    ids <- unique(upd$ID)
    after <- dplyr::tbl(con, "annotate") |>
      dplyr::filter(ID %in% !!ids) |>
      dplyr::select(ID, path, scaffold, annotate_lock) |>
      dplyr::collect()
    sel <- do.call(paste, after[UNIT_KEY]) %in% do.call(paste, upd[UNIT_KEY])
    after$annotate_lock[sel] <- 1L
    after <- after[after$annotate_lock %in% 1, , drop = FALSE]
    multi_path <- names(which(
      tapply(after$path, after$ID, function(p) length(unique(p))) > 1
    ))
    if (length(multi_path) > 0) {
      return(list(written = upd[0, , drop = FALSE], multi_path = multi_path))
    }
  }
  dplyr::tbl(con, "annotate") |>
    dplyr::rows_update(upd, unmatched = "ignore", in_place = TRUE, copy = TRUE, by = UNIT_KEY)
  list(written = upd, multi_path = character())
}

#' Is an export group name usable in export paths
#' @noRd
valid_group_name <- function(name) grepl("^[A-Za-z0-9._-]+$", name)

#' Assign units to an export group, or clear it with `NA`
#'
#' Upsert, not update: a unit that has never been grouped has no export row.
#'
#' @param group One name, a vector aligned to `units`, or `NA`.
#'
#' @return The data frame of rows upserted into `export`.
#' @noRd
assign_export_group <- function(con, units, group) {
  upd <- units[UNIT_KEY]
  stopifnot(length(group) %in% c(1L, nrow(upd)))
  upd$export_group <- rep_len(as.character(group), nrow(upd))
  dplyr::tbl(con, "export") |>
    dplyr::rows_upsert(upd, in_place = TRUE, copy = TRUE, by = UNIT_KEY)
  upd
}

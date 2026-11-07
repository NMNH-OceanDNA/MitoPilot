# R/user_annot_track.R
# Original (imported) user annotation track: place the user's features in the
# current sequence frame and flag what curation changed.

#' Reverse complement of an uppercase DNA string
#' @noRd
revcomp_chr <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
}

#' Locate `from` inside `to` by probe search
#'
#' @return list(offset, rc) such that position p of `from` lands at p + offset
#'   in `to` (or in reverseComplement(to) when rc), wrapped modulo nchar(to)
#'   when the lengths are equal. NULL when no probe hits exactly once.
#' @noRd
find_frame <- function(from, to, k = 60L) {
  Lf <- nchar(from)
  Lt <- nchar(to)
  if (Lf < k || Lt < k) return(NULL)
  same <- Lf == Lt
  targets <- list(to, revcomp_chr(to))
  if (same) targets <- lapply(targets, function(t) paste0(t, substr(t, 1L, k - 1L)))
  starts <- unique(pmax(1L, pmin(Lf - k + 1L, as.integer(c(1L, round(Lf * (1:7) / 8), Lf - k + 1L)))))
  for (p in starts) {
    probe <- substr(from, p, p + k - 1L)
    hits <- lapply(targets, function(t) {
      m <- gregexpr(probe, t, fixed = TRUE)[[1]]
      m[m > 0]
    })
    n <- lengths(hits)
    if (sum(n) == 0L) next
    if (sum(n) != 1L) next
    rc <- n[2] == 1L
    off <- as.integer(hits[[if (rc) 2L else 1L]][1]) - p
    if (same) {
      off <- off %% Lt
      if (off > 0) off <- off - Lt
    }
    return(list(offset = as.integer(off), rc = unname(rc)))
  }
  NULL
}

#' Map one feature through a frame
#'
#' Same length: circular shift (and mirror when rc). Different length (a trim):
#' linear shift, clipped to the sequence; an origin-spanning input is split at
#' the old origin and the larger surviving piece is kept.
#' @return list(pos1, pos2, dir, clipped) or NULL when nothing survives
#' @noRd
map_feature <- function(pos1, pos2, dir, frame, L_from, L_to) {
  flip <- function(d) if (identical(d, "-")) "+" else "-"
  mirror <- function(a, b) c(L_to - b + 1L, L_to - a + 1L)
  if (L_from == L_to) {
    q <- c(wrap_pos(pos1 + frame$offset, L_to), wrap_pos(pos2 + frame$offset, L_to))
    if (frame$rc) q <- mirror(q[1], q[2])
    return(list(pos1 = as.integer(q[1]), pos2 = as.integer(q[2]),
                dir = if (frame$rc) flip(dir) else dir, clipped = FALSE))
  }
  # A trim may follow a rotation (Linearize then Trim), and the probe gives one
  # offset for both, so the feature (unwrapped past the old origin when it
  # spans it) is tried shifted by -L, 0, +L and every copy that lands inside
  # the sequence is kept; the largest wins below.
  span <- if (pos1 > pos2) c(pos1, pos2 + L_from) else c(pos1, pos2)
  kept <- list()
  for (k in c(-2L, -1L, 0L, 1L)) {
    q <- span + frame$offset + k * L_from
    if (frame$rc) q <- mirror(q[1], q[2])
    a <- max(1L, q[1])
    b <- min(L_to, q[2])
    if (a > b) next
    kept[[length(kept) + 1L]] <- list(
      pos1 = as.integer(a), pos2 = as.integer(b), clipped = a != q[1] || b != q[2]
    )
  }
  if (!length(kept)) return(NULL)
  best <- kept[[which.max(vapply(kept, function(k) k$pos2 - k$pos1, numeric(1)))]]
  list(pos1 = best$pos1, pos2 = best$pos2,
       dir = if (frame$rc) flip(dir) else dir, clipped = best$clipped)
}

ORIG_NOTE_NO_SEQ <- "Original sequence not stored for this sample"
ORIG_NOTE_NO_PLACE <- "Original annotation cannot be placed on the edited sequence"

#' The sample's imported features placed in the unit's current frame
#'
#' Frame: original FASTA -> WF2 output (the assembly_backup sequence when an
#' app edit happened, else the current sequence) -> current sequence.
#' @return list(ok = TRUE, rows, length) or list(ok = FALSE, reason)
#' @noRd
original_features <- function(con, ID, path, scaffold) {
  fail <- function(r) list(ok = FALSE, reason = r)
  key <- list(as.character(ID), as.integer(path), as.integer(scaffold))
  if (!DBI::dbExistsTable(con, "user_contigs")) return(fail(ORIG_NOTE_NO_SEQ))
  uc <- DBI::dbGetQuery(con, "SELECT sequence FROM user_contigs WHERE ID = ?", params = key[1])
  if (nrow(uc) == 0 || is.na(uc$sequence[1])) return(fail(ORIG_NOTE_NO_SEQ))
  cur <- DBI::dbGetQuery(
    con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = key)
  if (nrow(cur) == 0 || is.na(cur$sequence[1])) return(fail(ORIG_NOTE_NO_PLACE))
  orig <- toupper(uc$sequence[1])
  cur <- toupper(cur$sequence[1])
  mid <- NULL
  if (has_assembly_backup(con, ID, path, scaffold)) {
    bk <- DBI::dbGetQuery(
      con, "SELECT sequence FROM assembly_backup WHERE ID = ? AND path = ? AND scaffold = ?",
      params = key)
    if (nrow(bk) == 0 || is.na(bk$sequence[1])) return(fail(ORIG_NOTE_NO_PLACE))
    mid <- toupper(bk$sequence[1])
  }
  wf2 <- mid %||% cur
  if (nchar(orig) != nchar(wf2)) return(fail(ORIG_NOTE_NO_PLACE))
  f1 <- find_frame(orig, wf2)
  if (is.null(f1)) return(fail(ORIG_NOTE_NO_PLACE))
  f2 <- if (!is.null(mid)) find_frame(mid, cur)
  if (!is.null(mid) && is.null(f2)) return(fail(ORIG_NOTE_NO_PLACE))

  ua <- DBI::dbGetQuery(con, "
    SELECT u.raw_name, u.feature_type, u.pos1, u.pos2, u.direction, u.segment_group, m.gene
    FROM user_annotations u LEFT JOIN gene_name_map m
      ON u.raw_name = m.raw_name AND u.feature_type = m.feature_type
    WHERE u.ID = ? ORDER BY u.feature_idx", params = key[1])
  L1 <- nchar(orig)
  L2 <- nchar(cur)
  rows <- lapply(seq_len(nrow(ua)), function(i) {
    m <- map_feature(ua$pos1[i], ua$pos2[i], ua$direction[i], f1, L1, L1)
    if (!is.null(f2)) m <- map_feature(m$pos1, m$pos2, m$dir, f2, L1, L2)
    if (is.null(m)) return(NULL)
    g <- ua$gene[i]
    type <- if (!is.na(g)) {
      if (g == "ctrl") "ctrl" else mito_gene_type(g)
    } else {
      switch(ua$feature_type[i], CDS = "PCG", tRNA = "tRNA", rRNA = "rRNA", "")
    }
    data.frame(
      name = ua$raw_name[i], gene = g, type = type,
      pos1 = m$pos1, pos2 = m$pos2, dir = m$dir,
      orig1 = as.integer(ua$pos1[i]), orig2 = as.integer(ua$pos2[i]), odir = ua$direction[i],
      clipped = m$clipped, segment_group = ua$segment_group[i],
      stringsAsFactors = FALSE
    )
  })
  list(ok = TRUE, rows = dplyr::bind_rows(rows), length = as.integer(L2))
}

#' Overlap in bp of two possibly origin-spanning intervals on a circle of length L
#' @noRd
circ_overlap_bp <- function(a1, a2, b1, b2, L) {
  pcs <- function(p1, p2) if (p1 > p2) list(c(p1, L), c(1L, p2)) else list(c(p1, p2))
  sum(vapply(pcs(a1, a2), function(x) {
    sum(vapply(pcs(b1, b2), function(y) max(0, min(x[2], y[2]) - max(x[1], y[1]) + 1), numeric(1)))
  }, numeric(1)))
}

#' Flag each placed original feature against the curated table
#'
#' Status, first that applies: unmapped, clipped, unchanged, moved, dropped.
#' Within a gene, originals and curated rows pair greedily by overlap bp, each
#' curated row claimed once.
#' @noRd
flag_original <- function(rows, curated, L) {
  n <- nrow(rows)
  status <- rep("dropped", n)
  if (n == 0) { rows$status <- character(0); return(rows) }
  status[is.na(rows$gene)] <- "unmapped"
  status[!is.na(rows$gene) & rows$clipped] <- "clipped"
  cur <- curated[!is.na(curated$pos1) & curated$pos1 > 0 & !is.na(curated$gene), , drop = FALSE]
  for (g in unique(rows$gene[status == "dropped"])) {
    oi <- which(rows$gene == g & status == "dropped")
    ci <- which(cur$gene == g)
    if (!length(ci)) next
    cand <- expand.grid(o = oi, c = ci)
    cand$bp <- mapply(function(o, c) circ_overlap_bp(rows$pos1[o], rows$pos2[o], cur$pos1[c], cur$pos2[c], L),
                      cand$o, cand$c)
    cand <- cand[cand$bp > 0, , drop = FALSE]
    cand <- cand[order(-cand$bp), , drop = FALSE]
    used_o <- integer(0); used_c <- integer(0)
    for (k in seq_len(nrow(cand))) {
      o <- cand$o[k]; c <- cand$c[k]
      if (o %in% used_o || c %in% used_c) next
      used_o <- c(used_o, o); used_c <- c(used_c, c)
      same <- rows$pos1[o] == cur$pos1[c] && rows$pos2[o] == cur$pos2[c] &&
        identical(rows$dir[o], cur$direction[c])
      status[o] <- if (same) "unchanged" else "moved"
    }
  }
  rows$status <- status
  rows
}

#' Original track rows as the viewer payload list, or NULL when there are none
#' @noRd
original_payload <- function(original) {
  if (!is.data.frame(original) || nrow(original) == 0) return(NULL)
  lapply(seq_len(nrow(original)), function(i) {
    f <- list(
      name = as.character(original$name[i]),
      gene = if (is.na(original$gene[i])) NULL else as.character(original$gene[i]),
      type = as.character(original$type[i]),
      pos1 = as.integer(original$pos1[i]), pos2 = as.integer(original$pos2[i]),
      dir = as.character(original$dir[i]), status = as.character(original$status[i]),
      orig1 = as.integer(original$orig1[i]), orig2 = as.integer(original$orig2[i]),
      odir = as.character(original$odir[i])
    )
    f[!vapply(f, is.null, logical(1))]
  })
}

#' One-line caption under a viewer's controls when the track cannot be shown
#' @noRd
orig_note_ui <- function(o) {
  if (is.null(o) || isTRUE(o$ok)) return(NULL)
  div(class = "mp-coverage-caption", o$reason)
}

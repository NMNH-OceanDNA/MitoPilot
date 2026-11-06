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
    if (sum(n) != 1L) return(NULL)
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
  pieces <- if (pos1 > pos2) list(c(pos1, L_from), c(1L, pos2)) else list(c(pos1, pos2))
  kept <- list()
  for (pc in pieces) {
    q <- pc + frame$offset
    if (frame$rc) q <- mirror(q[1], q[2])
    a <- max(1L, q[1])
    b <- min(L_to, q[2])
    if (a > b) next
    kept[[length(kept) + 1L]] <- list(
      pos1 = as.integer(a), pos2 = as.integer(b),
      clipped = a != q[1] || b != q[2] || length(pieces) > 1L
    )
  }
  if (!length(kept)) return(NULL)
  best <- kept[[which.max(vapply(kept, function(k) k$pos2 - k$pos1, numeric(1)))]]
  list(pos1 = best$pos1, pos2 = best$pos2,
       dir = if (frame$rc) flip(dir) else dir, clipped = best$clipped)
}

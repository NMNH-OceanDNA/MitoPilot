# R/app_annotate_synteny.R

#' Synteny viewer payload (tools/synteny_viewer_spec.md, section 4)
#'
#' Pure: no Shiny, no database. Coordinates pass through as stored; the browser
#' maps them to alignment columns. NULL when there is no usable alignment.
#' @noRd
synteny_payload <- function(annotations, blast_ref, aln, sample_label, ref_acc,
                            ref_label, selected = NULL, version = 0L) {
  if (is.null(aln) || nrow(aln) == 0) return(NULL)
  s <- toupper(as.character(aln$aligned_sample[1]))
  r <- toupper(as.character(aln$aligned_ref[1]))
  if (is.na(s) || is.na(r) || !nzchar(s) || !nzchar(r)) return(NULL)

  feat <- function(df, i, row = NULL) {
    f <- list(
      gene = as.character(df$gene[i]), type = as.character(df$type[i]),
      pos1 = as.integer(df$pos1[i]), pos2 = as.integer(df$pos2[i]),
      dir = as.character(df$direction[i])
    )
    if (!is.null(row)) f <- c(list(row = row), f)
    f
  }
  a <- annotations
  rows <- which(!is.na(a$pos1) & a$pos1 > 0 & !is.na(a$pos2))
  sample_feats <- lapply(rows, function(i) feat(a, i, row = i))
  ref_feats <- if (is.null(blast_ref) || nrow(blast_ref) == 0) list() else
    lapply(seq_len(nrow(blast_ref)), function(i) feat(blast_ref, i))
  sel <- if (length(selected) == 1L && !is.na(selected) && selected %in% rows) {
    as.integer(selected)
  }
  list(
    sampleLabel = sample_label, refAcc = ref_acc, refLabel = ref_label,
    refRc = isTRUE(aln$ref_rc[1]), s = s, r = r,
    refLength = as.integer(aln$ref_length[1]),
    rotation = as.integer(aln$rotation[1]),
    sampleFeatures = sample_feats, refFeatures = ref_feats,
    selected = sel, version = as.integer(version)
  )
}

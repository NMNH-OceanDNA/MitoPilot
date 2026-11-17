# R/app_annotate_synteny.R

#' Synteny viewer payload (tools/synteny_viewer_spec.md, section 4)
#'
#' Pure: no Shiny, no database. Coordinates pass through as stored; the browser
#' maps them to alignment columns. NULL when there is no usable alignment.
#' @noRd
synteny_payload <- function(annotations, blast_ref, aln, sample_label, ref_acc,
                            ref_label, selected = NULL, version = 0L,
                            extra_sample = NULL, hl_ref = NULL) {
  if (is.null(aln) || nrow(aln) == 0) return(NULL)
  s <- toupper(as.character(aln$aligned_sample[1]))
  r <- toupper(as.character(aln$aligned_ref[1]))
  if (is.na(s) || is.na(r) || !nzchar(s) || !nzchar(r)) return(NULL)

  feat <- function(df, i, row = NULL, hl = FALSE) {
    f <- list(
      gene = as.character(df$gene[i]), type = as.character(df$type[i]),
      pos1 = as.integer(df$pos1[i]), pos2 = as.integer(df$pos2[i]),
      dir = as.character(df$direction[i])
    )
    if (!is.null(row)) f <- c(list(row = row), f)
    if (hl) f$hl <- TRUE
    f
  }
  a <- annotations
  rows <- which(!is.na(a$pos1) & a$pos1 > 0 & !is.na(a$pos2))
  sample_feats <- lapply(rows, function(i) feat(a, i, row = i))
  if (!is.null(extra_sample)) {
    sample_feats <- c(sample_feats, lapply(seq_len(nrow(extra_sample)), function(i) {
      feat(extra_sample, i, row = as.character(extra_sample$row[i]), hl = TRUE)
    }))
  }
  ref_feats <- if (is.null(blast_ref) || nrow(blast_ref) == 0) list() else
    lapply(seq_len(nrow(blast_ref)), function(i) {
      feat(blast_ref, i, hl = !is.null(hl_ref) &&
             isTRUE(blast_ref$gene[i] == hl_ref$gene && blast_ref$pos1[i] == hl_ref$pos1))
    })
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

#' Synteny viewer controls + canvas block (tools/synteny_viewer_spec.md, section 5)
#' @noRd
synteny_viewer_ui <- function(ns, prefix = "synteny") {
  main <- prefix == "synteny"
  btn <- function(action, label, icon = NULL, title = NULL) {
    tags$button(
      type = "button", class = "btn btn-default", `data-mpsyn` = action,
      title = title, `aria-label` = title %||% label, icon, label
    )
  }
  tagList(
    div(
      class = "mp-seqview-controls",
      btn("fit", "Fit gene"),
      btn("whole", "Whole alignment"),
      btn("zoom_in", NULL, icon("magnifying-glass-plus"), "Zoom in"),
      btn("zoom_out", NULL, icon("magnifying-glass-minus"), "Zoom out"),
      if (main) uiOutput(ns("liftover_btn"), inline = TRUE)
    ),
    div(class = "mp-coverage-caption",
        if (main) {
          paste("Drag or scroll sideways to pan, scroll or pinch to zoom; click a sample gene to select its row,",
                "or a reference gene, then Lift over, to copy it to this assembly.",
                "Bases appear when zoomed in.")
        } else {
          "Drag or scroll sideways to pan, scroll or pinch to zoom. Bases appear when zoomed in."
        }),
    div(
      class = "mp-seqview",
      tags$canvas(id = ns(paste0(prefix, "_canvas")), class = "mp-seqview-canvas"),
      tags$div(id = ns(paste0(prefix, "_tip")), class = "mp-maptoref-tip", hidden = NA)
    )
  )
}

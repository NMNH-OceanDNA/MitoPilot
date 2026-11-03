#' Set partial_start / partial_stop from contig edges and codon evidence
#'
#' 5' (3') partial when the feature's 5' (3') end sits on the edge of a linear
#' contig, or for a PCG when the start codon (full 3-base stop codon) is not in
#' the gene rule. Rows sharing a JOIN group are one feature: only the 5'-most
#' member's start and the 3'-most member's stop are tested. Changed flags get
#' a note.
#'
#' @param annotations annotations data frame
#' @param assembly DNAStringSet
#' @param contig_key named character, contig id -> FASTA header
#' @param rules curation rules keyed by gene
#' @noRd
assess_partial_ends <- function(annotations, assembly, contig_key, rules) {
  n <- nrow(annotations)
  if (n == 0L) {
    return(annotations)
  }
  incoming <- function(col) {
    if (!col %in% names(annotations)) {
      return(integer(n))
    }
    dplyr::coalesce(as.integer(annotations[[col]]), 0L)
  }
  old_s <- incoming("partial_start")
  old_e <- incoming("partial_stop")

  hdr <- contig_key[annotations$contig]
  circ <- stringr::str_detect(hdr, "circular")
  len <- Biostrings::width(assembly[hdr])
  p1 <- annotations$pos1
  p2 <- annotations$pos2
  plus <- annotations$direction == "+"
  pcg <- annotations$type == "PCG"

  bad_start <- purrr::map2_lgl(annotations$gene, annotations$start_codon, \(g, cd) {
    opts <- rules[[g]][["start_codons"]]
    !is.null(opts) && !cd %in% opts
  })
  bad_stop <- purrr::map2_lgl(annotations$gene, annotations$stop_codon, \(g, cd) {
    opts <- rules[[g]][["stop_codons"]]
    !is.null(opts) && .codon_len(cd) == 3L && !cd %in% opts
  })

  new_s <- (!circ & ifelse(plus, p1 == 1L, p2 == len)) | (pcg & bad_start)
  new_e <- (!circ & ifelse(plus, p2 == len, p1 == 1L)) | (pcg & bad_stop)
  new_s <- as.integer(new_s)
  new_e <- as.integer(new_e)

  grp <- join_group_id(annotations$notes)
  for (g in unique(stats::na.omit(grp))) {
    m <- which(grp == g)
    # order members along the plus strand; on a circular contig the feature
    # starts after the largest gap between members
    ord <- m[order(p1[m])]
    if (circ[m[1]] && length(ord) > 1L) {
      gaps <- (p1[c(ord[-1], ord[1])] - p2[ord]) %% len[m[1]]
      first <- which.max(gaps) %% length(ord) + 1L
      ord <- ord[c(first:length(ord), seq_len(first - 1L))]
    }
    five <- if (plus[m[1]]) ord[1] else ord[length(ord)]
    three <- if (plus[m[1]]) ord[length(ord)] else ord[1]
    new_s[setdiff(m, five)] <- 0L
    new_e[setdiff(m, three)] <- 0L
  }

  note_at <- function(i, msg) {
    annotations$notes[i] <<- semicolon_paste(annotations$notes[i], msg)
  }
  note_at(which(new_s > old_s), "5' partial set by curation")
  note_at(which(new_s < old_s), "5' partial cleared by curation")
  note_at(which(new_e > old_e), "3' partial set by curation")
  note_at(which(new_e < old_e), "3' partial cleared by curation")

  annotations$partial_start <- new_s
  annotations$partial_stop <- new_e
  annotations
}

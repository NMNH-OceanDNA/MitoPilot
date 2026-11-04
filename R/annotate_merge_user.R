#' Drop tool annotation rows that collide with user-imported features
#'
#' @param tool tool annotation table (contig, gene, pos1, pos2)
#' @param user user feature table (contig, gene, pos1, pos2)
#' @param contig_lens named integer vector, contig -> length
#'
#' @noRd
#'
drop_overlapping_tool_rows <- function(tool, user, contig_lens) {
  if (nrow(tool) == 0L || nrow(user) == 0L) {
    return(tool)
  }
  ukey <- ifelse(user$gene == "OH", "ctrl", user$gene)
  tkey <- ifelse(tool$gene == "OH", "ctrl", tool$gene)
  drop <- purrr::pmap_lgl(
    list(tkey, tool$contig, tool$pos1, tool$pos2),
    \(g, ctg, p1, p2) {
      on_ctg <- user$contig == ctg
      hit <- on_ctg & circ_overlap(p1, p2, user$pos1, user$pos2)
      if (any(hit & ukey == g)) {
        return(TRUE)
      }
      hits <- user[hit, , drop = FALSE]
      if (nrow(hits) == 0L) {
        return(FALSE)
      }
      L <- contig_lens[[ctg]]
      total <- sum(purrr::map2_int(
        hits$pos1, hits$pos2, \(q1, q2) circ_overlap_len(p1, p2, q1, q2, L)
      ))
      total / circ_len(p1, p2, L) > 0.10
    }
  )
  tool[!drop, , drop = FALSE]
}

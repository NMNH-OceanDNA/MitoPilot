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
  key <- \(d) {
    d$gene <- ifelse(d$gene == "OH", "ctrl", d$gene)
    d
  }
  tool[!overlapped(key(tool), key(user), contig_lens, same_gene = TRUE), , drop = FALSE]
}

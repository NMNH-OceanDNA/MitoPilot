#' Add start/stop codons and translation for PCG rows
#'
#' Fills `start_codon`, `stop_codon`, `translation` for `type == "PCG"` rows
#' (NA elsewhere). Contigs match `names(assembly)` by first whitespace token.
#' Reads across the origin when pos1 > pos2. A `codon_start` column (1-3), if
#' present, shifts where translation begins on the coding strand.
#'
#' @param annotations tibble with contig, type, pos1, pos2, direction
#' @param assembly DNAStringSet
#' @param genetic_code NCBI genetic code id
#'
#' @noRd
add_pcg_codons <- function(annotations, assembly, genetic_code) {
  gc <- Biostrings::getGeneticCode(as.character(genetic_code))
  lookup <- stats::setNames(names(assembly), sub("\\s.*$", "", names(assembly)))
  added <- !"codon_start" %in% names(annotations)
  if (added) annotations$codon_start <- NA_integer_

  out <- annotations |>
    dplyr::rowwise() |>
    dplyr::mutate(
      .cds = if (type == "PCG" && !is.na(lookup[contig])) {
        suppressWarnings({
          s <- extract_circ_region(assembly[lookup[[contig]]], pos1, pos2)
          if (direction == "-") s <- Biostrings::reverseComplement(s)
          unname(as.character(s))
        })
      } else NA_character_,
      start_codon = if (!is.na(.cds)) substr(.cds, 1L, 3L) else NA_character_,
      stop_codon = if (!is.na(.cds)) {
        n <- nchar(.cds)
        len <- if (n %% 3L == 0L) 3L else n %% 3L
        substr(.cds, n - len + 1L, n)
      } else NA_character_,
      translation = if (!is.na(.cds)) {
        cs <- if (is.na(codon_start)) 1L else as.integer(codon_start)
        suppressWarnings(
          Biostrings::DNAString(substr(.cds, cs, nchar(.cds) - nchar(stop_codon))) |>
            Biostrings::translate(genetic.code = gc, if.fuzzy.codon = "solve") |>
            as.character()
        )
      } else NA_character_
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.cds)
  if (added) out$codon_start <- NULL
  out
}

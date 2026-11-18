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
  annotations$.cs <- rep_len(codon_offset(annotations$codon_start), nrow(annotations))

  annotations |>
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
        n <- nchar(.cds) - (.cs - 1L)
        len <- if (n %% 3L == 0L) 3L else n %% 3L
        substr(.cds, nchar(.cds) - len + 1L, nchar(.cds))
      } else NA_character_,
      translation = if (!is.na(.cds)) {
        suppressWarnings(
          Biostrings::DNAString(substr(.cds, .cs, nchar(.cds) - nchar(stop_codon))) |>
            Biostrings::translate(genetic.code = gc, if.fuzzy.codon = "solve") |>
            as.character()
        )
      } else NA_character_
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-.cds, -.cs)
}

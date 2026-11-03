#' Initialize a new user-annotation project database
#'
#' Builds the user-assembly schema with [new_db_userAsmb()] (mitogenome search,
#' circularization, and scaffold joining off), then parses each sample's
#' annotation file into the `user_annotations` table and the gene names found
#' into `gene_name_map`. The annotation tools are switched off in the default
#' `annotate_opts` set.
#'
#' @param db_path Path to the new database file
#' @param mapping_fn Path to the mapping CSV file. Must contain columns "ID",
#'   "Taxon", "Assembly", and "Annotation", plus "R1" and "R2" unless
#'   `no_raw_data = TRUE`.
#' @param mapping_id Column name of the mapping file to use as the primary key
#' @param assembly_path Directory holding the assembly FASTA files
#' @param annotation_path Directory holding the annotation files
#' @param no_raw_data (logical) Initialize a project with no raw reads (default =
#'   FALSE).
#' @param ... Additional arguments passed to [new_db_userAsmb()].
#' @return The import summary (invisibly).
#' @export
#'
new_db_userAnnot <- function(
    db_path = "./.sqlite",
    mapping_fn = NULL,
    mapping_id = "ID",
    assembly_path = NULL,
    annotation_path = NULL,
    no_raw_data = FALSE,
    ...) {
  mapping_fn <- mapping_fn %||% "./mapping.csv"
  mapping <- utils::read.csv(mapping_fn)
  mapping[[mapping_id]] <- as.character(mapping[[mapping_id]])
  .report_issues(
    check_user_annotations(mapping, assembly_path, annotation_path, .issues(),
                           labels = mapping[[mapping_id]]),
    "Annotation files"
  )

  new_db_userAsmb(
    db_path = db_path, mapping_fn = mapping_fn, mapping_id = mapping_id,
    assembly_path = assembly_path, no_raw_data = no_raw_data,
    find_mitogenome = FALSE, attempt_circularization = FALSE, join_scaffolds = FALSE,
    ...
  )

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(
    con,
    "CREATE TABLE user_annotations (
      ID TEXT NOT NULL, feature_idx INTEGER NOT NULL, seqid_raw TEXT, contig TEXT,
      feature_type TEXT, raw_name TEXT, raw_product TEXT,
      pos1 INTEGER, pos2 INTEGER, direction TEXT, anticodon TEXT,
      partial_start INTEGER DEFAULT 0, partial_stop INTEGER DEFAULT 0,
      codon_start INTEGER DEFAULT 1, segment_group INTEGER, segment_mode TEXT,
      segment_note TEXT, pseudo INTEGER DEFAULT 0, note TEXT, attrs TEXT,
      source_file TEXT, contig_length INTEGER,
      PRIMARY KEY (ID, feature_idx)
    );"
  )
  DBI::dbExecute(
    con,
    "CREATE TABLE gene_name_map (
      raw_name TEXT NOT NULL, feature_type TEXT NOT NULL DEFAULT '',
      gene TEXT, guess TEXT, manual INTEGER DEFAULT 0,
      n_samples INTEGER, n_features INTEGER, confirmed INTEGER DEFAULT 0,
      time_stamp TEXT,
      PRIMARY KEY (raw_name, feature_type)
    );"
  )
  DBI::dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_user_annotations_name
                       ON user_annotations(raw_name, feature_type)")
  imp <- import_user_annotations(con, mapping, assembly_path, annotation_path, mapping_id)
  DBI::dbExecute(
    con,
    "UPDATE annotate_opts SET use_mitos = 0, use_trnaScan = 0, use_arwen = 0,
       use_aragorn = 0, use_mitofinder = 0, coverage_trim = 0, start_gene = NULL
     WHERE annotate_opts = 'default'"
  )
  invisible(imp)
}

#' Parse every sample's annotation file into user_annotations and gene_name_map
#'
#' @return list(features = per-sample feature counts, dropped_types = per-sample
#'   named counts, warnings = per-sample parser warnings)
#' @noRd
import_user_annotations <- function(con, mapping, assembly_path, annotation_path,
                                    mapping_id = "ID") {
  ids <- as.character(mapping[[mapping_id]])
  parsed <- lapply(seq_along(ids), function(i) {
    L <- unname(Biostrings::fasta.seqlengths(.ua_path(assembly_path, mapping$Assembly[i])))[1]
    fn <- mapping$Annotation[i]
    p <- parse_user_annotation(.ua_path(annotation_path, fn), contig_length = L)
    ac <- sub("^([ACGTUacgtu]{3})$", "(\\1)", p$anticodon_raw)
    rows <- dplyr::tibble(
      ID = ids[i], feature_idx = seq_len(nrow(p)), seqid_raw = p$seqid,
      contig = paste0(ids[i], ".1.1"),
      feature_type = dplyr::coalesce(p$feature_type, ""),
      raw_name = dplyr::coalesce(p$raw_name, ""), raw_product = p$raw_product,
      pos1 = p$pos1, pos2 = p$pos2, direction = p$direction,
      anticodon = parse_anticodon(ac),
      partial_start = p$partial_start, partial_stop = p$partial_stop,
      codon_start = p$codon_start, segment_group = p$segment_group,
      segment_mode = p$segment_mode, segment_note = p$segment_note,
      pseudo = p$pseudo, note = p$note, attrs = p$attrs,
      source_file = basename(fn), contig_length = as.integer(L)
    )
    list(rows = rows, dropped = attr(p, "dropped_types"), warnings = attr(p, "warnings"))
  })
  ua <- dplyr::bind_rows(lapply(parsed, `[[`, "rows"))
  DBI::dbAppendTable(con, "user_annotations", ua)

  ua$guess <- guess_gene(ua$raw_name, dplyr::na_if(ua$feature_type, ""), ua$raw_product,
                         ua$note)$gene
  gm <- ua |>
    dplyr::group_by(raw_name, feature_type) |>
    dplyr::summarise(
      guess = guess[!is.na(guess)][1],
      n_samples = dplyr::n_distinct(ID),
      n_features = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(gene = guess, manual = 0L, confirmed = 0L, time_stamp = NA_character_)
  DBI::dbAppendTable(con, "gene_name_map", as.data.frame(gm))

  list(
    features = stats::setNames(vapply(parsed, function(x) nrow(x$rows), integer(1)), ids),
    dropped_types = stats::setNames(lapply(parsed, `[[`, "dropped"), ids),
    warnings = stats::setNames(lapply(parsed, `[[`, "warnings"), ids)
  )
}

#' Print imported feature counts and the gene-name table to the console
#'
#' @param import optional import summary from new_db_userAnnot(); adds dropped
#'   types and parser warnings, which the database does not keep
#' @noRd
user_annot_summary <- function(con, import = NULL) {
  ua <- DBI::dbGetQuery(con, "
    SELECT u.ID, u.feature_idx, u.segment_group, u.partial_start, u.partial_stop,
           u.pseudo, m.gene
    FROM user_annotations u LEFT JOIN gene_name_map m
      ON u.raw_name = m.raw_name AND u.feature_type = m.feature_type")
  ua$key <- ifelse(is.na(ua$segment_group), paste0("r", ua$feature_idx),
                   paste0("g", ua$segment_group))
  per <- ua |>
    dplyr::group_by(ID, key) |>
    dplyr::summarise(unmapped = all(is.na(gene)),
                     partial = any(partial_start == 1L | partial_stop == 1L),
                     pseudo = any(pseudo == 1L), .groups = "drop") |>
    dplyr::group_by(ID) |>
    dplyr::summarise(kept = dplyr::n(), unmapped = sum(unmapped),
                     partial = sum(partial), pseudo = sum(pseudo))
  lines <- "Imported annotations:"
  for (i in seq_len(nrow(per))) {
    id <- per$ID[i]
    dr <- import$dropped_types[[id]]
    drop_txt <- if (length(dr)) {
      paste0(sum(dr), " dropped (", paste(names(dr), dr, collapse = ", "), ")")
    } else if (!is.null(import)) "0 dropped"
    lines <- c(lines, paste0(
      "  ", id, ": ", paste(c(paste(per$kept[i], "kept"), drop_txt,
                              paste(per$unmapped[i], "unmapped"),
                              paste(per$partial[i], "partial"),
                              paste(per$pseudo[i], "pseudogene")), collapse = ", ")))
    for (w in import$warnings[[id]]) lines <- c(lines, paste0("    note: ", w))
  }
  gm <- DBI::dbGetQuery(con, "SELECT raw_name, feature_type, n_samples, n_features, guess
                              FROM gene_name_map ORDER BY guess IS NULL, guess, raw_name")
  gm$guess[is.na(gm$guess)] <- "(unmapped)"
  message(paste(c(lines, "", "Gene names found:",
                  utils::capture.output(print(gm, row.names = FALSE))), collapse = "\n"))
  invisible(gm)
}

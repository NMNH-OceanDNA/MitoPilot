# Package-default FASTA header templates for the Export modal. Used to seed the
# export_opts DB row and as the fallback when no custom template is stored.
# The leading token is {seqid}, not {ID}: a fragmented sample emits one record per
# scaffold and GenBank needs each SeqID unique within a submission. seqid is the
# plain ID for single-unit samples, so this is a no-op for unfragmented projects.
DEFAULT_FASTA_HEADER <- paste0(
  "{seqid} [organism={Taxon}] [topology={topology}] [mgcode={genetic_code}] ",
  "[location=mitochondrion] {Taxon} mitochondrion, {completeness}"
)
DEFAULT_FASTA_HEADER_GENE <- paste0(
  "{seqid} [organism={Taxon}] [topology={topology}] [mgcode={genetic_code}] ",
  "[location=mitochondrion] {Taxon}"
)
# The {ID}-era defaults, frozen for migration matching. Spelled out rather than
# derived from the current defaults: the default has changed twice ({completeness},
# then {seqid}), so deriving from the current value would stop matching the older
# ones. Only templates identical to a legacy default are migrated; custom templates
# are the user's own.
LEGACY_FASTA_HEADER_ID <- paste0(
  "{ID} [organism={Taxon}] [topology={topology}] [mgcode={genetic_code}] ",
  "[location=mitochondrion] {Taxon} mitochondrion, {completeness}"
)
LEGACY_FASTA_HEADER_GENE_ID <- paste0(
  "{ID} [organism={Taxon}] [topology={topology}] [mgcode={genetic_code}] ",
  "[location=mitochondrion] {Taxon}"
)

#' List saved export header template names
#'
#' Returns the names (the `export_opts` PK column) of saved templates, with
#' "default" guaranteed first. Falls back to just "default" when the table is
#' missing (pre-migration DB).
#'
#' @param con database connection
#'
#' @return character vector of template names
#'
#' @noRd
list_export_templates <- function(con) {
  names <- tryCatch(
    dplyr::tbl(con, "export_opts") |>
      dplyr::pull("export_opts"),
    error = function(e) character(0)
  )
  unique(c("default", sort(names)))
}

#' Read a saved FASTA header template from the project DB
#'
#' Returns the named template's strings, falling back to the package defaults
#' when the table/row is missing (e.g. a project DB that predates the
#' export_opts migration, or an unknown name).
#'
#' @param con database connection
#' @param name template name (defaults to "default")
#'
#' @return list(fasta_header=, fasta_header_gene=, start_aa=, stop_aa=,
#'   ident_pct=, aa_alignments=, gene_export=)
#'
#' @noRd
get_export_opts <- function(con, name = "default") {
  defaults <- list(
    fasta_header = DEFAULT_FASTA_HEADER,
    fasta_header_gene = DEFAULT_FASTA_HEADER_GENE,
    start_aa = 10, stop_aa = 10, ident_pct = 60, aa_alignments = TRUE, gene_export = FALSE
  )
  row <- tryCatch(
    dplyr::tbl(con, "export_opts") |>
      dplyr::filter(export_opts == !!name) |>
      dplyr::collect(),
    error = function(e) NULL
  )
  if (is.null(row) || nrow(row) == 0) return(defaults)
  pick <- function(x, default) {
    if (length(x) == 0 || is.na(x[1]) || !nzchar(x[1])) default else header_lock(x[1])
  }
  num <- function(k) {
    x <- row[[k]]
    if (length(x) == 0 || is.na(x[1])) defaults[[k]] else as.numeric(x[1])
  }
  list(
    fasta_header = pick(row$fasta_header, defaults$fasta_header),
    fasta_header_gene = pick(row$fasta_header_gene, defaults$fasta_header_gene),
    start_aa = num("start_aa"), stop_aa = num("stop_aa"), ident_pct = num("ident_pct"),
    aa_alignments = num("aa_alignments") == 1, gene_export = num("gene_export") == 1
  )
}

#' Save (upsert) a named FASTA header template and export options to the project DB
#'
#' Only the values given (not NULL or NA) are written. No-op (with a warning)
#' if the table is missing.
#'
#' @param con database connection
#' @param fasta_header,fasta_header_gene template strings
#' @param name template name (defaults to "default")
#' @param start_aa,stop_aa,ident_pct,aa_alignments,gene_export export options
#'
#' @noRd
set_export_opts <- function(con, fasta_header = NULL, fasta_header_gene = NULL, name = "default",
                            start_aa = NULL, stop_aa = NULL, ident_pct = NULL,
                            aa_alignments = NULL, gene_export = NULL) {
  vals <- list(fasta_header = fasta_header, fasta_header_gene = fasta_header_gene,
               start_aa = start_aa, stop_aa = stop_aa, ident_pct = ident_pct,
               aa_alignments = aa_alignments, gene_export = gene_export)
  vals <- Filter(function(v) length(v) == 1 && !is.na(v), vals)
  vals <- lapply(vals, function(v) if (is.logical(v)) as.integer(v) else v)
  tryCatch(
    dplyr::tbl(con, "export_opts") |>
      dplyr::rows_upsert(
        data.frame(export_opts = name, vals),
        in_place = TRUE,
        copy = TRUE,
        by = "export_opts"
      ),
    error = function(e) warning("Could not save export header templates: ", e$message)
  )
  invisible(NULL)
}

# Pre-export checks shared by the Export Data window and `mitopilot export`.
# `d` is the Export table data (export_table_data()).

export_table_data <- function(con) meta_view_join(fetch_export_data(con, session = NULL), con)

export_group_rows <- function(d, group) d[!is.na(d$export_group) & d$export_group == group, , drop = FALSE]

# Samples contributing units from more than one assembly path: blocked, since
# export_files() rejects them (each path would be a duplicate record).
export_multi_path_ids <- function(d, group) {
  d <- export_group_rows(d, group)
  if (nrow(d) == 0) return(character(0))
  by_id <- split(d$path, d$ID)
  names(by_id)[vapply(by_id, function(p) length(unique(p)) > 1, logical(1))]
}

# Samples whose one path is split across scaffolds, so each exports as several
# incomplete records: warned, not blocked.
export_fragmented_ids <- function(d, group) {
  d <- export_group_rows(d, group)
  if (nrow(d) == 0) return(character(0))
  by_id <- split(d, d$ID)
  names(by_id)[vapply(by_id, function(x) nrow(x) > 1 && length(unique(x$path)) == 1, logical(1))]
}

export_shown <- function(x, n) {
  s <- paste(utils::head(x, n), collapse = ", ")
  if (length(x) > n) paste0(s, ", and ", length(x) - n, " more") else s
}

export_multi_path_text <- function(mp) {
  shown <- export_shown(mp, 8)
  stringr::str_glue(
    "{mp_n(length(mp), 'sample')} still have more than one assembly path: {shown}.\n\n",
    "Assembly paths are alternative resolutions of the same genome, so ",
    "exporting each would submit duplicate records for one specimen. In the ",
    "Assemble module, open the assembly details and 'ignore' all but the ",
    "correct path (or build a consensus Path 0), then export again."
  )
}

export_fragmented_text <- function(frag) {
  shown <- export_shown(frag, 5)
  stringr::str_glue(
    "{mp_n(length(frag), 'sample')} have more than one assembly and will each ",
    "produce a SEPARATE GenBank record: {shown}.\n\n",
    "That is correct when the scaffolds really are different genomes. If a ",
    "sample is instead ONE genome broken into fragments, each record will ",
    "be submitted as an incomplete genome. Cancel and use consensus ",
    "trimming / scaffold joining to combine them, or 'ignore' all but one ",
    "scaffold."
  )
}

export_nmnh_text <- function(r) {
  sprintf(paste(
    "%d of %s in this group have NMNH problems. NMNH records should not be",
    "submitted to GenBank without both a specimen_voucher and a voucherURI.",
    "Exporting anyway leaves the empty modifiers out of the headers."),
    sum(!r$ok), mp_n(nrow(r), "sample"))
}

# Export data with the NMNH values from nmnh_resolve() (NA when not resolved)
export_with_nmnh <- function(d, r) {
  for (k in c("nmnh_specimen_voucher", "nmnh_voucherURI")) {
    d[[k]] <- if (is.null(r)) rep(NA_character_, nrow(d)) else r[[k]][match(d$ID, r$ID)]
  }
  d
}

# Message for the first invalid header template, or NULL. The gene header is
# only checked when genes are exported.
export_header_problem <- function(main, gene, gene_export, d) {
  v_main <- validate_fasta_header(main, d)
  v_gene <- if (isTRUE(gene_export)) validate_fasta_header(gene, d) else list(ok = TRUE)
  if (isTRUE(v_main$ok) && isTRUE(v_gene$ok)) return(NULL)
  bad <- if (!isTRUE(v_main$ok)) v_main$message else v_gene$message
  which_t <- if (!isTRUE(v_main$ok)) "mitogenome" else "gene"
  stringr::str_glue("The {which_t} header template is invalid: {bad}")
}

# Conflicts on the specimen items the header templates use. "not checked" items
# never count, and a template without specimen tokens never warns.
export_specimen_rows <- function(con, ids, templates) {
  concepts <- specimen_template_concepts(templates, specimen_csv_columns(con))
  if (!length(concepts)) return(NULL)
  specimen_export_warnings(specimen_conflicts(con, ids), concepts, ids)
}

# A short "...near here..." context window around a brace position, for messages.
brace_snippet <- function(chars, pos) {
  n <- length(chars)
  lo <- max(1L, pos - 6L)
  hi <- min(n, pos + 10L)
  paste0(
    if (lo > 1L) "..." else "",
    paste(chars[lo:hi], collapse = ""),
    if (hi < n) "..." else ""
  )
}

#' Find the first unmatched brace in a template
#'
#' Treats every open and close brace literally (no glue escape special-casing)
#' and requires them to balance: each close brace matches the most recent open
#' brace. Reports the first brace that has no match: a close brace with no open
#' brace before it, or, if all close, the first open brace left at the end.
#'
#' @param template header template string
#'
#' @return NULL when braces balance, else
#'   list(pos = integer, type = "open"|"close", snippet = character)
#'
#' @noRd
find_unmatched_brace <- function(template) {
  if (length(template) != 1 || is.na(template)) return(NULL)
  chars <- strsplit(template, "", fixed = TRUE)[[1]]
  n <- length(chars)
  open_stack <- integer(0)  # positions of '{' awaiting a match
  for (i in seq_len(n)) {
    c <- chars[i]
    if (c == "{") {
      open_stack <- c(open_stack, i)
    } else if (c == "}") {
      if (length(open_stack) == 0L) {
        return(list(pos = i, type = "close", snippet = brace_snippet(chars, i)))
      }
      open_stack <- open_stack[-length(open_stack)]
    }
  }
  if (length(open_stack) > 0L) {
    pos <- open_stack[1]  # first '{' left unclosed
    return(list(pos = pos, type = "open", snippet = brace_snippet(chars, pos)))
  }
  NULL
}

#' Validate a glue-syntax FASTA header template
#'
#' Checks brace balance, then dry-runs the template through
#' [stringr::str_glue_data()] against a representative data row so malformed
#' braces and unknown `{column}` references are caught before they reach
#' [export_files()] (where they would crash the GUI). Returns a cleaned,
#' user-facing message rather than the raw glue error.
#'
#' @param template header template string
#' @param data data frame whose columns the template may reference (e.g. rv$data).
#'   When empty, the template is still parsed for brace balance.
#' @param require_completeness when TRUE, an otherwise-valid template that does
#'   not end with `{completeness}` returns a non-blocking warning (ok = TRUE,
#'   level = "warn") so the user can still save/export.
#'
#' @return list(ok = logical, level = "ok"|"warn"|"error", message = character).
#'   Blocking errors (unbalanced braces, unknown column, empty) return
#'   `ok = FALSE`.
#'
#' @noRd
validate_fasta_header <- function(template, data = NULL, require_completeness = FALSE) {
  err <- function(msg) list(ok = FALSE, level = "error", message = msg)
  if (is.null(template) || !nzchar(trimws(template))) {
    return(err("Template is empty."))
  }
  template <- header_clean(template)
  # Braces must balance: every { needs a matching }. Point at the first
  # unmatched one so the user can find it.
  brace <- find_unmatched_brace(template)
  if (!is.null(brace)) {
    b <- if (brace$type == "open") "{" else "}"
    return(err(sprintf(
      "unmatched '%s' at position %d (near \"%s\")",
      b, brace$pos, brace$snippet
    )))
  }
  row <- if (!is.null(data) && nrow(data) > 0) data[1, , drop = FALSE] else data.frame()
  tok <- regmatches(template, gregexpr("\\{[^{}]*\\}", template))[[1]]
  unk <- setdiff(substr(tok, 2L, nchar(tok) - 1L), names(row))
  if (ncol(row) && length(unk)) {
    return(err(sprintf(
      '{%s} is not one of the available fields. Pick one from the "Available metadata" panel.', unk[1]
    )))
  }
  tryCatch({
    stringr::str_glue_data(row, template)
    has_data <- !is.null(data) && nrow(data) > 0
    hdr <- if (has_data) header_clean(stringr::str_glue_data(data, template)) else character()
    fmt <- header_format_check(template, hdr, data$ID)
    if (identical(fmt$level, "error")) return(fmt)
    # A template without {seqid} gives every unit of a multi-assembly sample the
    # same defline, while the .tbl still carries the per-unit >Feature seqid. The
    # two must agree exactly or table2asn rejects the submission. seqid is the
    # plain ID for single-unit samples, so {seqid} is always the safe choice.
    if (!grepl("\\{seqid\\}", template)) {
      multi_unit <- !is.null(data) && "ID" %in% names(data) && any(duplicated(data$ID))
      msg <- paste(
        "Header does not use {seqid}. Samples with more than one assembly",
        "will produce duplicate FASTA deflines that do not match the .tbl",
        ">Feature line, and table2asn will reject the submission. Use {seqid}",
        "instead of {ID}: it is the plain ID for samples with one assembly."
      )
      if (multi_unit) {
        return(err(msg))
      }
      return(list(ok = TRUE, level = "warn", message = msg))
    }
    if (has_data) {
      dm <- unique(unlist(duplicate_modifiers(hdr)))
      nm <- dm[gsub("_", "", dm) %in% c("specimenvoucher", "voucheruri")]
      if (length(nm)) {
        return(err(sprintf("Duplicate %s. NMNH records need each exactly once: take the other out.",
                           paste0("[", nm, "=]", collapse = ", "))))
      }
      if (length(dm)) {
        return(list(ok = TRUE, level = "warn", message = sprintf(
          paste("Duplicate %s. NCBI accepts each modifier once: keep one",
                "and take the other out."),
          paste0("[", dm, "=]", collapse = ", "))))
      }
    }
    if (!is.null(fmt)) return(fmt)
    if (has_data) {
      mf <- missing_fields(template, data)
      if (!is.null(mf)) return(list(ok = TRUE, level = "warn", message = mf))
    }
    if (require_completeness && !grepl("\\{completeness\\}\\s*$", template)) {
      return(list(
        ok = TRUE, level = "warn",
        message = paste(
          "Header does not end with {completeness}; GenBank submissions may not",
          "reflect partial vs complete genome status. You can export anyway."
        )
      ))
    }
    list(ok = TRUE, level = "ok", message = "Valid template.")
  }, error = function(e) {
    raw <- conditionMessage(e)
    # Unknown field -> glue reports "object 'XXX' not found"
    col <- regmatches(raw, regexpr("object '[^']+' not found", raw))
    if (length(col) > 0) {
      name <- sub("object '([^']+)' not found", "\\1", col)
      return(err(sprintf(
        '{%s} is not one of the available fields. Pick one from the "Available metadata" panel.',
        name
      )))
    }
    # Fallback: strip glue's multi-line wrapper to the last informative line
    err(sub("^.*!\\s*", "", gsub("\n", " ", raw)))
  })
}

#' Message naming header columns that some samples leave empty, or NULL
#' @noRd
missing_fields <- function(template, data) {
  tok <- regmatches(template, gregexpr("\\{[^{}]*\\}", template))[[1]]
  nm <- intersect(unique(substr(tok, 2L, nchar(tok) - 1L)), names(data))
  nm <- nm[!startsWith(nm, "nmnh_")]
  ids <- unique(data$ID)
  out <- vapply(nm, function(n) {
    v <- trimws(as.character(data[[n]]))
    miss <- unique(data$ID[is.na(v) | !nzchar(v) | v == "NA"])
    if (!length(miss)) return(NA_character_)
    sprintf("{%s} for %d of %d sample(s) (%s)", n, length(miss), length(ids),
            paste(utils::head(miss, 3), collapse = ", "))
  }, character(1))
  out <- out[!is.na(out)]
  if (length(out)) {
    paste0("Missing data in ", paste(out, collapse = "; "),
           ". A [modifier=] with no value is left out of that sample's header.")
  }
}

#' Names of the source modifiers used more than once in each header
#' @noRd
duplicate_modifiers <- function(headers) {
  m <- regmatches(headers, gregexpr("\\[[^]=[]+=", headers))
  lapply(m, function(x) {
    nm <- gsub("[- ]", "_", tolower(trimws(substr(x, 2L, nchar(x) - 1L))))
    unique(nm[duplicated(nm)])
  })
}

#' GenBank SeqID for an assembly unit
#'
#' Plain `ID` when the sample contributes exactly one exported record, otherwise
#' `ID_p<path>_s<scaffold>`. `n_units` counts only exported (non-ignored) units, so
#' the suffix appears only where there is a sibling record to disambiguate from and
#' single-scaffold projects keep their existing names.
#'
#' @param ID,path,scaffold unit key (vectorised).
#' @param n_units number of exported units for that unit's sample.
#'
#' @noRd
export_seqid <- function(ID, path, scaffold, n_units) {
  # Recycle n_units explicitly: ifelse() returns a result the length of its test, so
  # a scalar n_units (e.g. dplyr::n() inside a grouped mutate) would otherwise
  # collapse every unit to the first one's SeqID.
  n <- rep_len(as.integer(n_units), length(ID))
  ifelse(
    n > 1L,
    paste0(ID, "_p", path, "_s", scaffold),
    as.character(ID)
  )
}

#' Populate export table
#'
#' @param db database connection
#' @param session reactive session
#' @param locked_only keep only units whose assembly and annotation are locked
#'   (the Export tab's view). `FALSE` returns every non-ignored unit, for callers
#'   that describe units they have already exported.
#'
#' @noRd
fetch_export_data <- function(con = NULL, session = getDefaultReactiveDomain(),
                              locked_only = TRUE) {
  db <- con %||% session$userData$con

  samples <- dplyr::tbl(db, "samples") |>
    dplyr::select(-dplyr::any_of("topology"))

  # ORF count per annotate unit, blanked when ORF finding is disabled. annotate and
  # annotations are keyed (ID, path, scaffold), so both the grouping and the joins
  # below must carry the full key: keying on ID alone would sum ORFs across a
  # sample's scaffolds and fan each row out to one copy per scaffold.
  unit_key <- c("ID", "path", "scaffold")
  orf_counts <- dplyr::tbl(db, "annotations") |>
    dplyr::select(ID, path, scaffold, type) |>
    dplyr::collect() |>
    dplyr::group_by(ID, path, scaffold) |>
    dplyr::summarise(ORFCount = sum(type == "ORF", na.rm = TRUE), .groups = "drop")
  orf_enabled <- dplyr::tbl(db, "annotate") |>
    dplyr::select(ID, path, scaffold, orf_opts) |>
    dplyr::left_join(dplyr::tbl(db, "orf_opts"), by = "orf_opts") |>
    dplyr::select(ID, path, scaffold, use_orffinder) |>
    dplyr::collect()

  # Resolved per-unit reference (override, else this scaffold's own BLAST top hit),
  # the same source the Annotate table, synteny default and export note use. Each
  # unit is independent: a scaffold with no hit must not inherit a sibling's
  # accession, which is what reading the sample-level assemble.blast_* did.
  assemblies_unit <- unit_ref_facts(db) |>
    dplyr::select(ID, path, scaffold, blast_accession, blast_accession_auto,
                  blast_species, blast_lineage)

  # Topology is per scaffold. assemblies carries this contig's own value, while
  # annotate.topology can summarise a whole sample ("fragmented", or a legacy
  # joined string), and only a per-scaffold value belongs in the Export table or
  # the summary CSV. Matches what export_files() writes to the defline.
  unit_topology <- dplyr::tbl(db, "assemblies") |>
    dplyr::filter(ignore != 1) |>
    dplyr::select(ID, path, scaffold, scaffold_topology = topology) |>
    dplyr::collect()

  out <- dplyr::tbl(db, "assemble")
  if (locked_only) out <- dplyr::filter(out, assemble_lock == 1)
  out <- out |>
    dplyr::select(ID, dplyr::any_of("poor_blast_ref")) |>
    dplyr::left_join(dplyr::tbl(db, "annotate"), by = "ID")
  if (locked_only) out <- dplyr::filter(out, annotate_lock == 1)
  out <- out |>
    dplyr::select(
      ID, path, scaffold, curate_opts, topology,
      length, structure, PCGCount, tRNACount, rRNACount, missing, extra, warnings,
      dplyr::any_of(c("annotate_switch", "blast_accession_auto",
                      "poor_blast_ref", "partial"))
    ) |>
    dplyr::left_join(
      dplyr::tbl(db, "curate_opts") |>
        dplyr::select(curate_opts, dplyr::any_of("linear_complete")),
      by = "curate_opts"
    ) |>
    dplyr::left_join(samples, by = "ID") |>
    # Export state is per unit; a unit with no row yet has never been grouped.
    dplyr::left_join(dplyr::tbl(db, "export"), by = unit_key) |>
    dplyr::select(-R1, -R2) |>
    dplyr::relocate(Taxon, .after = ID) |>
    dplyr::collect() |>
    .specimen_status_join(db) |>
    dplyr::relocate(specimen, .after = Taxon) |>
    # inner_join: gates rows to a non-ignored assembly, as assemblies_unit is
    # already filtered on ignore.
    dplyr::inner_join(assemblies_unit, by = unit_key) |>
    dplyr::relocate(blast_accession, blast_accession_auto, blast_species,
                    blast_lineage, .after = specimen) |>
    dplyr::left_join(orf_counts, by = unit_key) |>
    dplyr::left_join(orf_enabled, by = unit_key) |>
    dplyr::left_join(unit_topology, by = unit_key)

  out <- .meta_join(out, db)

  # these columns are absent on un-migrated DBs
  if (!"linear_complete" %in% names(out)) out$linear_complete <- NA_integer_
  if (!"partial" %in% names(out)) out$partial <- NA_character_
  if (!"annotate_switch" %in% names(out)) out$annotate_switch <- NA_integer_

  out |>
    dplyr::mutate(
      # Per-scaffold topology wins, and anything not usable as a topology is
      # coerced to "linear". Same rule and same order as export_files(), so the
      # table can never disagree with the defline it is previewing.
      topology = dplyr::if_else(
        !is.na(scaffold_topology) & scaffold_topology != "",
        scaffold_topology, topology
      ),
      topology = dplyr::if_else(
        !is.na(topology) & topology %in% c("circular", "linear"),
        topology, "linear"
      ),
      # Ref-align status is only meaningful when this unit has a real BLAST hit
      # (mirrors the Annotate table).
      blast_ref_status = dplyr::if_else(
        is.na(blast_accession) | blast_accession == "NO HIT",
        NA_character_, poor_blast_ref
      ),
      # Auto-derive completeness from topology; per-sample "partial" forces
      # partial; project-level linear_complete forces linear -> complete.
      completeness = dplyr::case_when(
        !is.na(partial) & partial == "yes" ~ "partial genome",
        topology == "circular" ~ "complete genome",
        !is.na(linear_complete) & linear_complete == 1L ~ "complete genome",
        TRUE ~ "partial genome"
      ),
      structure = stringr::str_replace_all(structure, "trn[A-Z]", "\u2022"),
      export_group = as.character(export_group),
      ORFCount = dplyr::if_else(
        is.na(use_orffinder) | use_orffinder != 1L,
        NA_integer_,
        as.integer(ORFCount)
      )
    ) |>
    dplyr::select(-use_orffinder, -scaffold_topology,
                  -dplyr::any_of("linear_complete")) |>
    # The SeqID this unit will export under. Counted within (ID, export_group) to
    # mirror export_files(): a sample is only suffixed when it contributes more than
    # one record to the same submission. Also what {seqid} resolves to when the
    # Export modal previews a header template.
    dplyr::group_by(ID, export_group) |>
    dplyr::mutate(seqid = export_seqid(ID, path, scaffold, dplyr::n())) |>
    dplyr::ungroup() |>
    dplyr::relocate(path, scaffold, seqid, .after = ID) |>
    dplyr::relocate(ORFCount, .after = rRNACount) |>
    dplyr::relocate(blast_ref_status, .after = blast_accession) |>
    # State leads, like the other tables; Export Group is pinned right, so
    # nothing may render after it (T08).
    dplyr::relocate(annotate_switch, .before = ID) |>
    dplyr::relocate(export_group, .after = dplyr::last_col())
}

#' Mapping-file columns to show in the Export table
#'
#' Everything the user put in the mapping file that MitoPilot does not own or
#' already display. fetch_export_data() carries every `samples` column, so the
#' `samples` schema alone decides the set. Read paths and the user-assembly
#' inputs are MitoPilot's.
#'
#' @param sample_cols Column names of the `samples` table.
#' @param declared Column names the table already defines.
#' @return Character vector, in `samples` order.
#' @noRd
export_metadata_cols <- function(sample_cols, declared) {
  owned <- c("ID", "Taxon", "genetic_code", "topology", "R1", "R2", "assembly",
             "Assembly", "Annotation", "Topology", "Reference", "Reference_topology",
             "specimen", "specimen_message", "specimen_icons")
  sample_cols[!sample_cols %in% c(owned, declared)]
}

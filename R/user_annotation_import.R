#' Parse a user-supplied annotation file (GFF3, BED, or 5-column feature table)
#'
#' One row per segment. `pos1 > pos2` marks a feature spanning the origin.
#' Attributes: `dropped_types` (type -> count), `is_circular`, `warnings`.
#'
#' @param fn path to a .gff/.gff3, .bed, or .tbl file
#' @param contig_length contig length, needed to recognize origin-wrapping parts
#'
#' @noRd
parse_user_annotation <- function(fn, contig_length = NULL) {
  if (grepl("\\.gff3?$", fn, ignore.case = TRUE)) return(.read_gff3(fn, contig_length))
  if (grepl("\\.bed$", fn, ignore.case = TRUE)) return(.read_bed(fn))
  if (grepl("\\.tbl$", fn, ignore.case = TRUE)) return(.read_tbl(fn, contig_length))
  stop("Unrecognized annotation file type (expected .gff, .gff3, .bed, or .tbl): ", fn,
       call. = FALSE)
}

#' Value of a GFF3 attribute, case-insensitive key, URL-decoded
#'
#' @param attrs raw column 9 strings
#' @param key attribute name
#'
#' @noRd
.gff_attr <- function(attrs, key) {
  vapply(attrs, function(s) .ua_get(.gff_attr_list(s), key), character(1),
         USE.NAMES = FALSE)
}

.gff_attr_list <- function(s) {
  if (is.na(s) || s %in% c("", ".")) return(character(0))
  kv <- trimws(strsplit(s, ";", fixed = TRUE)[[1]])
  kv <- kv[nzchar(kv)]
  val <- ifelse(grepl("=", kv, fixed = TRUE), sub("^[^=]*=", "", kv), "")
  dec <- function(x) vapply(gsub("%(?![0-9A-Fa-f]{2})", "%25", x, perl = TRUE),
                            utils::URLdecode, character(1), USE.NAMES = FALSE)
  stats::setNames(dec(val), dec(sub("=.*$", "", kv)))
}

# first present of `keys` in a named character vector, case-insensitive
.ua_get <- function(a, keys) {
  for (k in keys) {
    i <- match(tolower(k), tolower(names(a)))
    if (!is.na(i)) return(unname(a[[i]]))
  }
  NA_character_
}

.ua_has <- function(a, key) tolower(key) %in% tolower(names(a))

.ua_drop_keys <- function(a, keys) a[!tolower(names(a)) %in% tolower(keys)]

.ua_lines <- function(fn) {
  txt <- sub("\\s+$", "", readLines(fn, warn = FALSE))
  data.frame(line = seq_along(txt), text = txt, stringsAsFactors = FALSE)
}

.ua_split <- function(text, n) {
  p <- lapply(strsplit(text, "\t", fixed = TRUE), function(x) {
    x <- trimws(x)
    length(x) <- n
    x
  })
  matrix(as.character(unlist(p)), ncol = n, byrow = TRUE)
}

.ua_one_seqid <- function(seqid, fn) {
  u <- unique(seqid)
  if (length(u) > 1L) {
    stop("Annotation file has more than one seqid (", paste(u, collapse = ", "),
         "); one contig per file: ", fn, call. = FALSE)
  }
  u
}

.ua_count <- function(x) {
  tb <- table(x)
  stats::setNames(as.integer(tb), as.character(names(tb)))
}

.ua_json <- function(v) {
  if (!length(v)) return("{}")
  g <- split(unname(v), factor(names(v), unique(names(v))))
  as.character(jsonlite::toJSON(lapply(g, unique), auto_unbox = TRUE))
}

.ua_clean_note <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) NA_character_ else gsub("[\t\r\n]", " ", gsub(";", ",", paste(x, collapse = ", ")))
}

.ua_note_name <- function(x) {
  x <- trimws(x)
  long <- !is.na(x) & nchar(x) > 40L
  cut <- substr(x[long], 1L, 41L)
  y <- sub("[[:space:][:punct:]]*\\s\\S*$", "", cut)
  y[y == cut] <- substr(cut[y == cut], 1L, 40L)
  x[long] <- paste0(y, "...")
  x
}

.ua_feat <- function(type, raw_name, segs, direction = "+", raw_product = NA_character_,
                     anticodon_raw = NA_character_, p5 = FALSE, p3 = FALSE,
                     cs_attr = NA_integer_, slip = FALSE, exception = NA_character_,
                     note = NA_character_, pseudo = FALSE, attrs = character(0),
                     feature_id = NA_character_, source_row = NA_integer_) {
  as.list(environment())
}

.ua_empty <- function() {
  dplyr::tibble(
    seqid = character(), feature_type = character(), raw_name = character(),
    raw_product = character(), pos1 = integer(), pos2 = integer(),
    direction = character(), anticodon_raw = character(),
    partial_start = integer(), partial_stop = integer(), codon_start = integer(),
    segment_group = integer(), segment_mode = character(), segment_note = character(),
    pseudo = integer(), note = character(), attrs = character(),
    feature_id = character(), source_row = integer()
  )
}

# shared by all readers: wrap merge, JOIN groups, 5'/3' flags, transl_except check
.ua_finish <- function(feats, seqid, contig_length, dropped, is_circular) {
  L <- contig_length
  grp <- 0L
  warn <- character(0)
  rows <- lapply(feats, function(f) {
    s <- f$segs
    minus <- identical(f$direction, "-")
    if (!is.null(L) && nrow(s) > 1L) {
      e <- which(s$pos2 == L)[1]
      o <- which(s$pos1 == 1L)[1]
      if (!is.na(e) && !is.na(o) && e != o) {
        if (minus) s$phase[e] <- s$phase[o]
        s$pos2[e] <- s$pos2[o]
        s <- s[-o, , drop = FALSE]
      }
    }
    key <- s$pos1
    w <- s$pos1 > s$pos2
    if (any(w)) key <- ifelse(!w & s$pos1 <= max(s$pos2[w]), key + L, key)
    i5 <- if (minus) which.max(key) else which.min(key)
    i3 <- if (minus) which.min(key) else which.max(key)
    ph <- suppressWarnings(as.integer(s$phase[i5]))
    cs <- if (!is.na(ph) && ph %in% 0:2) ph + 1L else
      if (!is.na(f$cs_attr)) as.integer(f$cs_attr) else 1L
    n <- nrow(s)
    idx <- seq_len(n)
    g <- NA_integer_
    mode <- NA_character_
    snote <- NA_character_
    if (n > 1L) {
      grp <<- grp + 1L
      g <- grp
      mode <- if (f$slip) "frameshift" else "exon"
      exc <- f$exception
      if (!is.na(exc) && grepl("^ribosomal slippage$", exc, ignore.case = TRUE)) exc <- NA
      snote <- .ua_clean_note(c(exc, f$note))
    }
    te <- .ua_get(f$attrs, "transl_except")
    len <- sum(ifelse(w, L - s$pos1 + 1L + s$pos2, s$pos2 - s$pos1 + 1L))
    if (!is.na(te) && grepl("aa:TERM", te, ignore.case = TRUE) && len %% 3 == 0) {
      warn <<- c(warn, paste0(
        "transl_except present but coordinates give a full stop codon: ", f$raw_name))
    }
    ac <- f$anticodon_raw
    if (is.na(ac)) {
      m <- stringr::str_match(f$raw_name, "\\(([acgtuACGTU]{3})\\)$|-([ACGTU]{3})$")
      ac <- dplyr::coalesce(m[, 2], m[, 3])
    }
    dplyr::tibble(
      seqid = seqid, feature_type = as.character(f$type), raw_name = f$raw_name,
      raw_product = f$raw_product, pos1 = as.integer(s$pos1), pos2 = as.integer(s$pos2),
      direction = f$direction, anticodon_raw = ac,
      partial_start = as.integer(idx == i5 & (f$p5 || cs > 1L)),
      partial_stop = as.integer(idx == i3 & f$p3),
      codon_start = ifelse(idx == i5, cs, 1L),
      segment_group = g, segment_mode = mode, segment_note = snote,
      pseudo = as.integer(f$pseudo), note = f$note, attrs = .ua_json(f$attrs),
      feature_id = f$feature_id, source_row = as.integer(f$source_row)
    )
  })
  out <- dplyr::bind_rows(.ua_empty(), unname(rows))
  structure(out, dropped_types = dropped, is_circular = is_circular, warnings = warn)
}

.ua_gff_keep <- c("cds", "trna", "rrna", "gene", "mrna", "exon", "trna_gene",
                  "rrna_gene", "ncrna_gene", "d_loop", "d-loop", "misc_feature",
                  "rep_origin", "origin_of_replication", "repeat_region", "pseudogene")

.ua_gff_consumed <- c("ID", "Parent", "Name", "gene", "product", "anticodon", "Note",
                      "partial", "start_range", "end_range", "codon_start", "pseudo",
                      "pseudogene", "exception", "ribosomal_slippage")

#' @noRd
.read_gff3 <- function(fn, contig_length = NULL) {
  ln <- .ua_lines(fn)
  fa <- grep("^##FASTA", ln$text)
  if (length(fa)) ln <- ln[seq_len(fa[1] - 1L), , drop = FALSE]
  ln <- ln[nzchar(ln$text) & !startsWith(ln$text, "#"), , drop = FALSE]
  cols <- .ua_split(ln$text, 9L)
  seqid <- .ua_one_seqid(cols[, 1], fn)
  lt <- tolower(cols[, 3])
  a <- lapply(cols[, 9], .gff_attr_list)

  reg <- which(lt == "region")
  ic <- if (length(reg)) .ua_get(a[[reg[1]]], "Is_circular") else NA_character_
  is_circular <- if (is.na(ic)) NA else tolower(ic) == "true"

  keep <- lt %in% .ua_gff_keep
  dropped_raw <- cols[!keep & lt != "region", 3]
  k <- which(keep)
  cols <- cols[k, , drop = FALSE]
  lt <- lt[k]
  a <- a[k]
  line <- ln$line[k]
  start <- as.integer(cols[, 4])
  end <- as.integer(cols[, 5])

  id <- vapply(a, .ua_get, character(1), "ID")
  par <- sub(",.*$", "", vapply(a, .ua_get, character(1), "Parent"))
  up <- match(par, id)
  is_parent <- lt %in% c("gene", "mrna", "pseudogene") | grepl("_gene$", lt)
  is_child <- lt %in% c("cds", "trna", "rrna")
  span <- paste(start, end)
  span_par <- which(is_parent)[match(span, span[is_parent])]
  up_child <- ifelse(is_child & is.na(up), span_par, up)

  gene <- vapply(a, .ua_get, character(1), "gene")
  pseudo <- lt == "pseudogene" | vapply(a, function(x) .ua_has(x, "pseudo") ||
                                          .ua_has(x, "pseudogene"), logical(1))
  has_child <- logical(length(lt))
  for (i in which(is_child)) {
    j <- up_child[i]
    seen <- integer(0)
    while (!is.na(j) && !j %in% seen) {
      seen <- c(seen, j)
      has_child[j] <- TRUE
      if (is.na(gene[i]) && !is.na(gene[j])) gene[i] <- gene[j]
      if (is.na(gene[i])) gene[i] <- .ua_get(a[[j]], "Name")
      if (lt[j] == "pseudogene") pseudo[i] <- TRUE
      j <- up[j]
    }
  }
  exon_abs <- lt == "exon" & !is.na(par)
  dropped <- .ua_count(c(dropped_raw, cols[exon_abs, 3]))
  r <- which(!(is_parent & has_child) & !exon_abs)

  gkey <- ifelse(!is.na(id), paste0("I:", id),
                 ifelse(lt == "cds" & !is.na(par), paste0("P:", par), paste0("R:", seq_along(lt))))[r]
  feats <- lapply(split(r, factor(gkey, unique(gkey))), function(ix) {
    i <- ix[1]
    ga <- a[ix]
    get1 <- function(key) {
      v <- vapply(ga, .ua_get, character(1), key)
      v[!is.na(v)][1]
    }
    nm <- dplyr::coalesce(gene[i], get1("Name"), get1("product"), id[i],
                          .ua_note_name(get1("Note")))
    txt <- c(vapply(ga, .ua_get, character(1), "Name"), vapply(ga, .ua_get, character(1), "Note"))
    minus <- cols[i, 7] %in% "-"
    lo_open <- any(vapply(ga, .ua_has, logical(1), "start_range"))
    hi_open <- any(vapply(ga, .ua_has, logical(1), "end_range"))
    exc <- get1("exception")
    extra <- c(source = cols[i, 2], score = cols[i, 6])
    extra <- extra[!is.na(extra) & extra != "."]
    .ua_feat(
      type = if (lt[i] == "pseudogene") "CDS" else cols[i, 3],
      raw_name = nm,
      segs = data.frame(pos1 = start[ix], phase = cols[ix, 8], pos2 = if (is.null(contig_length))
        end[ix] else ifelse(end[ix] > contig_length, end[ix] - contig_length, end[ix])),
      direction = if (minus) "-" else "+",
      raw_product = get1("product"),
      anticodon_raw = get1("anticodon"),
      p5 = (if (minus) hi_open else lo_open) ||
        any(grepl("5'\\s*partial", txt, ignore.case = TRUE)),
      p3 = (if (minus) lo_open else hi_open) ||
        any(grepl("3'\\s*partial", txt, ignore.case = TRUE)),
      cs_attr = suppressWarnings(as.integer(get1("codon_start"))),
      slip = any(vapply(ga, .ua_has, logical(1), "ribosomal_slippage")) ||
        (!is.na(exc) && grepl("ribosomal slippage|RNA editing", exc, ignore.case = TRUE)),
      exception = exc,
      note = .ua_clean_note(vapply(ga, .ua_get, character(1), "Note")),
      pseudo = any(pseudo[ix]),
      attrs = c(extra, unlist(lapply(ga, .ua_drop_keys, .ua_gff_consumed))),
      feature_id = id[i],
      source_row = line[i]
    )
  })
  .ua_finish(feats, seqid, contig_length, dropped, is_circular)
}

#' @noRd
.read_bed <- function(fn) {
  ln <- .ua_lines(fn)
  ln <- ln[nzchar(ln$text) & !grepl("^(#|track|browser)", ln$text), , drop = FALSE]
  cols <- .ua_split(ln$text, 12L)
  seqid <- .ua_one_seqid(cols[, 1], fn)
  if (any(is.na(cols[, 4]) | !nzchar(cols[, 4]))) {
    stop("BED file needs a name in column 4 on every row: ", fn, call. = FALSE)
  }
  feats <- lapply(seq_len(nrow(cols)), function(i) {
    r <- cols[i, ]
    s0 <- as.integer(r[2])
    nb <- suppressWarnings(as.integer(r[10]))
    segs <- if (!is.na(nb) && nb > 1L) {
      sz <- as.integer(strsplit(r[11], ",", fixed = TRUE)[[1]][seq_len(nb)])
      st <- as.integer(strsplit(r[12], ",", fixed = TRUE)[[1]][seq_len(nb)])
      data.frame(pos1 = s0 + st + 1L, pos2 = s0 + st + sz, phase = NA)
    } else {
      data.frame(pos1 = s0 + 1L, pos2 = as.integer(r[3]), phase = NA)
    }
    ex <- c(score = r[5], thickStart = r[7], thickEnd = r[8], itemRgb = r[9])
    .ua_feat(type = NA_character_, raw_name = r[4], segs = segs,
             direction = if (r[6] %in% "-") "-" else "+",
             attrs = ex[!is.na(ex)], feature_id = r[4], source_row = ln$line[i])
  })
  .ua_finish(feats, seqid, NULL, .ua_count(character(0)), NA)
}

.ua_tbl_keep <- c("gene", "CDS", "tRNA", "rRNA", "D-loop", "misc_feature", "rep_origin",
                  "origin_of_replication", "repeat_region")

.ua_tbl_consumed <- c("gene", "product", "note", "anticodon", "codon_start", "pseudo",
                      "pseudogene", "exception", "ribosomal_slippage")

#' @noRd
.read_tbl <- function(fn, contig_length = NULL) {
  ln <- .ua_lines(fn)
  seqid <- NA_character_
  ft <- list()
  for (i in seq_len(nrow(ln))) {
    s <- ln$text[i]
    n <- length(ft)
    if (grepl("^>Feature\\s", s)) {
      if (!is.na(seqid)) stop("Feature table has more than one seqid: ", fn, call. = FALSE)
      seqid <- sub("\\|$", "", sub("^[a-z]+\\|", "", sub("^>Feature\\s+(\\S+).*$", "\\1", s)))
    } else if (grepl("^[<>]?\\d+\t[<>]?\\d+\t\\S", s)) {
      p <- strsplit(s, "\t", fixed = TRUE)[[1]]
      ft[[n + 1L]] <- list(key = trimws(p[3]), iv = list(p[1:2]), q = character(0),
                           line = ln$line[i])
    } else if (n && grepl("^[<>]?\\d+\t[<>]?\\d+$", s)) {
      ft[[n]]$iv <- c(ft[[n]]$iv, list(strsplit(s, "\t", fixed = TRUE)[[1]]))
    } else if (n && grepl("^\t\t\t\\S", s)) {
      p <- strsplit(sub("^\t\t\t", "", s), "\t", fixed = TRUE)[[1]]
      v <- if (length(p) > 1L) trimws(paste(p[-1], collapse = "\t")) else ""
      ft[[n]]$q <- c(ft[[n]]$q, stats::setNames(v, trimws(p[1])))
    }
  }

  num <- function(x) as.integer(gsub("[<>]", "", x))
  ends <- vapply(ft, function(f) paste(num(f$iv[[1]][1]), num(f$iv[[length(f$iv)]][2])),
                 character(1))
  mark <- function(g, x) if (grepl("^[<>]", g) && !grepl("^[<>]", x)) paste0(substr(g, 1, 1), x) else x
  keys <- vapply(ft, `[[`, character(1), "key")
  drop <- !keys %in% .ua_tbl_keep
  for (i in which(keys == "gene" & seq_along(ft) < length(ft))) {
    j <- i + 1L
    if (keys[j] %in% c("CDS", "tRNA", "rRNA") && ends[i] == ends[j]) {
      drop[i] <- TRUE
      gi <- ft[[i]]$iv
      nj <- length(ft[[j]]$iv)
      ft[[j]]$iv[[1]][1] <- mark(gi[[1]][1], ft[[j]]$iv[[1]][1])
      ft[[j]]$iv[[nj]][2] <- mark(gi[[length(gi)]][2], ft[[j]]$iv[[nj]][2])
      if (!.ua_has(ft[[j]]$q, "gene") && .ua_has(ft[[i]]$q, "gene")) {
        ft[[j]]$q <- c(ft[[j]]$q, gene = .ua_get(ft[[i]]$q, "gene"))
      }
    }
  }
  dropped <- .ua_count(keys[!keys %in% .ua_tbl_keep])

  feats <- lapply(ft[!drop], function(f) {
    q <- f$q
    iv <- f$iv
    d <- Filter(function(v) v[1] != v[2], lapply(iv, num))
    notes <- unname(q[names(q) == "note"])
    ac_note <- grepl("^anticodon:", notes)
    p5_note <- grepl("5' end not determined", notes, ignore.case = TRUE)
    p3_note <- grepl("3' end not determined", notes, ignore.case = TRUE)
    rest <- notes[!(ac_note | p5_note | p3_note)]
    exc <- .ua_get(q, "exception")
    .ua_feat(
      type = f$key,
      raw_name = dplyr::coalesce(.ua_get(q, "gene"), .ua_get(q, "product"),
                                 .ua_note_name(rest[1])),
      segs = data.frame(
        pos1 = vapply(iv, function(v) min(num(v)), integer(1)),
        pos2 = vapply(iv, function(v) max(num(v)), integer(1)),
        phase = NA
      ),
      direction = if (length(d) && d[[1]][1] > d[[1]][2]) "-" else "+",
      raw_product = .ua_get(q, "product"),
      anticodon_raw = dplyr::coalesce(.ua_get(q, "anticodon"), notes[ac_note][1]),
      p5 = grepl("^[<>]", iv[[1]][1]) || any(p5_note),
      p3 = grepl("^[<>]", iv[[length(iv)]][2]) || any(p3_note),
      cs_attr = suppressWarnings(as.integer(.ua_get(q, "codon_start"))),
      slip = .ua_has(q, "ribosomal_slippage") ||
        (!is.na(exc) && grepl("ribosomal slippage|RNA editing", exc, ignore.case = TRUE)),
      exception = exc,
      note = .ua_clean_note(rest),
      pseudo = .ua_has(q, "pseudo") || .ua_has(q, "pseudogene"),
      attrs = .ua_drop_keys(q, .ua_tbl_consumed),
      source_row = f$line
    )
  })
  .ua_finish(feats, seqid, contig_length, dropped, NA)
}

#' Map a GFF/tbl feature type to a MitoPilot gene type (NA if ambiguous).
#' @noRd
.gff_type_to_mito <- function(feature_type) {
  ft <- as.character(feature_type)
  out <- rep(NA_character_, length(ft))
  out[ft %in% c("CDS", "mRNA")] <- "PCG"
  out[ft %in% c("tRNA", "tRNA_gene")] <- "tRNA"
  out[ft %in% c("rRNA", "rRNA_gene")] <- "rRNA"
  out[tolower(ft) %in% c("d_loop", "d-loop", "rep_origin", "origin_of_replication")] <- "ctrl"
  out
}

#' Guess canonical gene name and type for imported features. Vectorised.
#' @noRd
guess_gene <- function(raw_name, feature_type, raw_product = NA, raw_note = NA) {
  n <- max(length(raw_name), length(feature_type), length(raw_product), length(raw_note))
  raw_name <- rep_len(as.character(raw_name), n)
  raw_product <- rep_len(as.character(raw_product), n)
  raw_note <- rep_len(as.character(raw_note), n)
  mt <- rep_len(.gff_type_to_mito(feature_type), n)
  hit <- function(p, x = list(raw_name, raw_product, raw_note)) {
    Reduce(`|`, lapply(x, grepl, pattern = p, ignore.case = TRUE, perl = TRUE))
  }
  cr_re <- paste0("control region|d-loop|d loop|displacement loop|",
                  "^(D[-_ ]?loop|CR|OH|control[ _]region|rep_origin)$")
  cr <- mt %in% "ctrl"
  ol <- hit("L-strand|light[- ]strand|\\bOL\\b|origin of L")
  mt[cr] <- NA_character_
  one <- function(nm, pr, ty) {
    tys <- if (is.na(ty)) c("PCG", "rRNA", "tRNA") else ty
    for (t in tys) {
      g <- normalize_mito_gene(nm, t, pr)
      if (!is.na(g) && !startsWith(g, "?")) return(g)
    }
    NA_character_
  }
  gene <- mapply(function(nm, pr, ty) {
    g <- one(nm, pr, ty)
    tok <- sub("[[:space:]_].*$", "", nm)
    if (is.na(g) && !is.na(tok) && nzchar(tok) && tok != nm &&
        !grepl("orf|intron|like", nm, ignore.case = TRUE)) g <- one(tok, pr, ty)
    g
  }, raw_name, raw_product, mt, USE.NAMES = FALSE)
  gene <- as.character(gene)
  keep <- !is.na(gene) & gene %in% MITO_GENE_CHOICES
  gene[!keep | cr] <- NA_character_
  ctrl <- !ol & (cr | hit(cr_re, list(raw_name)) |
                   (is.na(gene) & hit(cr_re, list(raw_product, raw_note))))
  gene[ctrl] <- "ctrl"
  type <- vapply(gene, function(g) if (is.na(g)) NA_character_ else if (g == "ctrl") "ctrl" else
    mito_gene_type(g), character(1), USE.NAMES = FALSE)
  dplyr::tibble(gene = gene, type = type)
}

#' Extract an uppercase DNA anticodon from NCBI, export-note, or name forms.
#' @noRd
parse_anticodon <- function(x) {
  one <- function(s) {
    if (is.na(s)) return(NA_character_)
    pats <- c("seq:([acgtuACGTU]{3})", "anticodon:([acgtuACGTU]{3})",
              "\\(([acgtuACGTU]{3})\\)\\s*$", "-([ACGTU]{3})$")
    for (p in pats) {
      m <- regmatches(s, regexec(p, s))[[1]]
      if (length(m) == 2) return(chartr("U", "T", toupper(m[2])))
    }
    NA_character_
  }
  vapply(as.character(x), one, character(1), USE.NAMES = FALSE)
}

#' Read an imported annotation CSV as annotate()'s merged table
#'
#' Keeps rows on contigs present in `assembly`, checks coordinates against
#' contig widths, and fills PCG codons and translation.
#'
#' @param fn path to the CSV
#' @param assembly DNAStringSet
#' @param genetic_code NCBI genetic code id
#'
#' @noRd
.read_user_annotations <- function(fn, assembly, genetic_code) {
  ann <- readr::read_csv(fn, col_types = readr::cols(.default = "c"), na = c("", "NA")) |>
    dplyr::mutate(dplyr::across(
      dplyr::any_of(c("pos1", "pos2", "length", "partial_start", "partial_stop", "codon_start")),
      as.integer
    ))
  lens <- stats::setNames(Biostrings::width(assembly), sub("\\s.*$", "", names(assembly)))
  ann <- ann[ann$contig %in% names(lens), ]
  if (nrow(ann) == 0L) {
    stop("user annotations: no features for contig ", paste(names(lens), collapse = ", "),
         " in ", fn, call. = FALSE)
  }
  for (ctg in unique(ann$contig)) {
    rows <- ann$contig == ctg
    n <- sum(pmax(ann$pos1[rows], ann$pos2[rows]) > lens[[ctg]], na.rm = TRUE)
    if (n > 0) {
      stop("user annotations: ", n, " feature(s) extend past the end of contig ", ctg,
           call. = FALSE)
    }
  }
  add_pcg_codons(ann, assembly, genetic_code) |>
    dplyr::arrange(contig, pos1)
}

.ua_trn_aa <- function(gene) unname(trnA_key_MITOS[sub("^(trn[LS])$", "\\11", gene)])

#' Imported features for one sample as annotate()'s user-annotation CSV
#'
#' Attribute `warnings`: JOIN groups whose members map to different genes.
#' @noRd
resolve_user_annotations <- function(con, ID) {
  ua <- DBI::dbGetQuery(con, "
    SELECT u.*, m.gene FROM user_annotations u LEFT JOIN gene_name_map m
      ON u.raw_name = m.raw_name AND u.feature_type = m.feature_type
    WHERE u.ID = ? ORDER BY u.feature_idx", params = list(ID))
  ua <- ua |>
    dplyr::group_by(segment_group) |>
    dplyr::mutate(split = !is.na(segment_group) & dplyr::n_distinct(gene) > 1L) |>
    dplyr::ungroup()
  bad <- ua[ua$split, ]
  warn <- vapply(split(bad, bad$segment_group), function(b) paste0(
    ID, ": segments of ", paste(unique(b$raw_name), collapse = "/"), " map to different genes (",
    paste(unique(dplyr::coalesce(b$gene, "dropped")), collapse = ", "),
    "); kept as separate features"), character(1), USE.NAMES = FALSE)
  ua <- ua[!is.na(ua$gene), ]
  g <- ifelse(ua$split, NA, ua$segment_group)
  g <- match(g, unique(g[!is.na(g)]))
  join <- ifelse(is.na(g), NA, paste0(
    "JOIN: mode=", ua$segment_mode, " group=", g,
    ifelse(is.na(ua$segment_note), "", paste0(" note=", ua$segment_note))))
  imp <- ifelse(!is.na(ua$note) & nzchar(ua$note) &
                  !(!is.na(join) & dplyr::coalesce(ua$note == ua$segment_note, FALSE)),
                paste0("import: ", ua$note), NA)
  parts <- cbind(join, ifelse(ua$pseudo %in% 1L, "import: pseudogene", NA), imp)
  notes <- apply(parts, 1, function(x) Reduce(semicolon_paste, x[!is.na(x)], NA_character_))
  type <- vapply(ua$gene, function(g) if (g == "ctrl") "ctrl" else mito_gene_type(g),
                 character(1), USE.NAMES = FALSE)
  aa <- .ua_trn_aa(ua$gene)
  out <- dplyr::tibble(
    contig = ua$contig, type = type, gene = ua$gene,
    product = dplyr::case_when(
      type == "PCG" ~ unname(CDS_key[ua$gene]),
      ua$gene == "rrnL" ~ "16S ribosomal RNA",
      ua$gene == "rrnS" ~ "12S ribosomal RNA",
      type == "tRNA" & !is.na(aa) ~ paste0("tRNA-", aa),
      type == "ctrl" ~ "d-loop",
      .default = NA_character_
    ),
    pos1 = ua$pos1, pos2 = ua$pos2,
    length = ifelse(ua$pos1 > ua$pos2, ua$contig_length - ua$pos1 + 1L + ua$pos2,
                    ua$pos2 - ua$pos1 + 1L),
    direction = ua$direction, anticodon = ua$anticodon,
    partial_start = ua$partial_start, partial_stop = ua$partial_stop,
    codon_start = ua$codon_start, notes = as.character(notes), tool = rep("user", nrow(ua))
  ) |>
    dplyr::arrange(contig, pos1)
  structure(out, warnings = warn)
}

#' Write every sample's resolved CSV; returns IDs whose file content changed
#'
#' Attribute `warnings`: collected from resolve_user_annotations().
#' @noRd
write_user_annotation_csvs <- function(con, project_dir) {
  ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples")$ID
  dir <- file.path(project_dir, "user_annotations")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  warn <- character(0)
  changed <- vapply(ids, function(id) {
    r <- resolve_user_annotations(con, id)
    warn <<- c(warn, attr(r, "warnings"))
    txt <- readr::format_csv(r, na = "")
    fn <- file.path(dir, paste0(id, "_annotations.csv"))
    old <- if (file.exists(fn)) readChar(fn, file.size(fn), useBytes = TRUE) else NA_character_
    if (identical(old, txt)) return(FALSE)
    writeChar(txt, fn, eos = NULL)
    TRUE
  }, logical(1), USE.NAMES = FALSE)
  structure(ids[changed], warnings = warn)
}

#' TRUE when every gene_name_map row is confirmed
#' @noRd
user_annot_confirmed <- function(con) {
  x <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n, SUM(confirmed = 1) AS k FROM gene_name_map")
  x$n > 0 && identical(as.numeric(x$n), as.numeric(x$k))
}

#' Confirm the gene-name mapping, write CSVs, requeue changed samples
#'
#' @return list(kept, dropped, requeued, warnings)
#' @noRd
confirm_user_annotations <- function(con, project_dir) {
  changed <- write_user_annotation_csvs(con, project_dir)
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "UPDATE gene_name_map SET confirmed = 1, time_stamp = ?",
                   params = list(as.numeric(Sys.time())))
    if (length(changed)) {
      DBI::dbExecute(con, paste0(
        "UPDATE annotate SET annotate_switch = 1 WHERE annotate_lock = 0 AND ID IN (",
        paste(rep("?", length(changed)), collapse = ", "), ")"), params = as.list(changed))
    }
  })
  n <- DBI::dbGetQuery(con, "
    SELECT SUM(m.gene IS NOT NULL) AS kept, SUM(m.gene IS NULL) AS dropped
    FROM user_annotations u LEFT JOIN gene_name_map m
      ON u.raw_name = m.raw_name AND u.feature_type = m.feature_type")
  list(kept = as.integer(n$kept), dropped = as.integer(n$dropped),
       requeued = as.character(changed), warnings = attr(changed, "warnings"))
}

#' Set the gene for a raw name; NA or "(drop)" drops it
#' @noRd
set_gene_name_map <- function(con, raw_name, feature_type, gene) {
  gene <- as.character(gene)
  gene[is.na(gene) | gene == "(drop)"] <- NA_character_
  bad <- !is.na(gene) & !gene %in% MITO_IMPORT_GENE_CHOICES
  if (any(bad)) stop("set_gene_name_map: unknown gene '", gene[bad][1], "'")
  DBI::dbExecute(con, "UPDATE gene_name_map SET gene = ?, manual = 1, confirmed = 0
                       WHERE raw_name = ? AND feature_type = ?",
                 params = list(gene, raw_name, dplyr::coalesce(feature_type, "")))
}

#' Features carrying one imported name, one row per feature
#' @noRd
gene_name_map_samples <- function(con, raw_name, feature_type) {
  DBI::dbGetQuery(con, "
    SELECT ID, feature_type, COUNT(*) OVER (PARTITION BY ID) AS n_features, pos1, pos2,
           direction, partial_start, partial_stop, codon_start, segment_group, pseudo,
           note, source_file
    FROM user_annotations WHERE raw_name = ? AND feature_type = ?
    ORDER BY ID, pos1", params = list(raw_name, dplyr::coalesce(feature_type, "")))
}

#' Restore every gene to its automatic guess
#' @noRd
reset_gene_name_map <- function(con) {
  DBI::dbExecute(con, "UPDATE gene_name_map SET gene = guess, manual = 0, confirmed = 0")
}

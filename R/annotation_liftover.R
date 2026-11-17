#' Liftover thresholds (defaults and hard floors)
#' @noRd
LIFTOVER_DEFAULTS <- list(min_pctid = 85, max_gap_pct = 5)
LIFTOVER_FLOORS <- list(min_pctid = 80, max_gap_pct = 10)
LIFTOVER_PICK_MSG <- "Select a reference gene in the plot to lift over to the sample assembly."

liftover_clean <- function(s) gsub("[^ACGTNacgtnRYSWKMBDHVryswkmbdhv]", "N", s)

#' Gene-orientation sequence of a reference feature (spliced when pos1 > pos2)
#' @noRd
liftover_ref_gene_seq <- function(ref_seq, pos1, pos2, direction) {
  s <- Biostrings::DNAString(liftover_clean(ref_seq))
  L <- length(s)
  as.character(feature_nt(s, wrap_pos(pos1, L), wrap_pos(pos2, L), direction))
}

#' Reason a reference feature cannot be lifted over, or NA
#' @noRd
liftover_blocked_reason <- function(blast_ref, gene, pos1) {
  type <- blast_ref$type[blast_ref$gene == gene & blast_ref$pos1 == pos1][1]
  if (identical(type, "ctrl")) return("Control regions cannot be lifted over.")
  if (identical(type, "PCG") && sum(blast_ref$gene == gene & blast_ref$type == "PCG", na.rm = TRUE) > 1) {
    return("Split gene models cannot be lifted over.")
  }
  NA_character_
}

#' Map a reference gene onto an assembly; hits passing LIFTOVER_FLOORS, best first
#' @noRd
liftover_map <- function(gene_seq, assembly_seq, circular, max_hits = 3) {
  out <- data.frame(
    pos1 = integer(), pos2 = integer(), direction = character(),
    pctid = numeric(), gap_pct = numeric(), score = numeric(),
    aln_gene = character(), aln_sample = character()
  )
  gene_seq <- liftover_clean(gene_seq)
  assembly_seq <- liftover_clean(assembly_seq)
  n <- nchar(gene_seq)
  L <- nchar(assembly_seq)
  if (!circular && n > L) return(out)
  g <- Biostrings::DNAString(gene_seq)
  g_rc <- Biostrings::reverseComplement(g)
  subj <- Biostrings::DNAString(
    if (circular) paste0(assembly_seq, substr(assembly_seq, 1, n)) else assembly_seq
  )
  subMx <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -1, baseOnly = FALSE)
  fit <- function(pat) {
    pwalign::pairwiseAlignment(
      pattern = pat, subject = subj, substitutionMatrix = subMx,
      gapOpening = 10, gapExtension = 4, type = "global-local"
    )
  }
  count_gaps <- function(x) nchar(gsub("[^-]", "", as.character(x)))
  for (k in seq_len(max_hits)) {
    fwd <- fit(g)
    rev <- fit(g_rc)
    minus <- pwalign::score(rev) > pwalign::score(fwd)
    a <- if (minus) rev else fwd
    s <- BiocGenerics::start(pwalign::subject(a))
    e <- BiocGenerics::end(pwalign::subject(a))
    ap <- pwalign::alignedPattern(a)
    asub <- pwalign::alignedSubject(a)
    x <- strsplit(as.character(ap), "")[[1]]
    y <- strsplit(as.character(asub), "")[[1]]
    pctid <- 100 * sum(x == y & x %in% c("A", "C", "G", "T")) / length(x)
    gap_pct <- 100 * (count_gaps(ap) + count_gaps(asub)) / n
    if (pctid < LIFTOVER_FLOORS$min_pctid || gap_pct > LIFTOVER_FLOORS$max_gap_pct) break
    if (minus) {
      ap <- Biostrings::reverseComplement(Biostrings::DNAString(as.character(ap)))
      asub <- Biostrings::reverseComplement(Biostrings::DNAString(as.character(asub)))
    }
    out[k, ] <- list(
      wrap_pos(s, L), wrap_pos(e, L), if (minus) "-" else "+",
      pctid, gap_pct, pwalign::score(a), as.character(ap), as.character(asub)
    )
    f <- wrap_pos(s:e, L)
    idx <- unique(c(f, f + L))
    idx <- idx[idx <= length(subj)]
    subj <- Biostrings::replaceLetterAt(subj, idx, strrep("N", length(idx)))
  }
  out
}

#' Shared bases between a hit and each live annotation (strand-blind, wrap-aware)
#' @noRd
liftover_overlaps <- function(pos1, pos2, annotations, asm_len) {
  ann <- annotations[!grepl("_DELETED_", annotations$gene, fixed = TRUE), , drop = FALSE]
  i <- seq_len(nrow(ann))
  bp <- vapply(i, function(j) circ_overlap_len(pos1, pos2, ann$pos1[j], ann$pos2[j], asm_len), numeric(1))
  ann_len <- vapply(i, function(j) circ_len(ann$pos1[j], ann$pos2[j], asm_len), numeric(1))
  hit_len <- circ_len(pos1, pos2, asm_len)
  out <- data.frame(
    gene = ann$gene, bp = bp,
    substantial = bp > 0.25 * hit_len | bp > 0.25 * ann_len
  )
  out[out$bp > 0, , drop = FALSE]
}

#' Clamp user thresholds to LIFTOVER_FLOORS
#' @noRd
liftover_clamp <- function(min_pctid, max_gap_pct) {
  list(
    min_pctid = max(min_pctid, LIFTOVER_FLOORS$min_pctid),
    max_gap_pct = min(max_gap_pct, LIFTOVER_FLOORS$max_gap_pct)
  )
}

#' Product name for a lifted gene, mirroring the annotation tools
#' @noRd
liftover_product <- function(gene, type, template) {
  p <- template$product[template$gene == gene & !is.na(template$product)][1]
  if (!is.na(p)) return(p)
  if (gene == "rrnL") return("16S ribosomal RNA")
  if (gene == "rrnS") return("12S ribosomal RNA")
  if (type == "tRNA" && !is.na(trnA_key_MITOS[gene])) return(paste0("tRNA-", trnA_key_MITOS[[gene]]))
  if (type == "PCG") return(unname(CDS_key[gene]))
  NA_character_
}

LIFTOVER_ANN_COLS <- c(
  "ID", "path", "scaffold", "type", "gene", "product", "pos1", "pos2", "length",
  "direction", "anticodon", "tool", "start_codon", "stop_codon", "partial_start",
  "partial_stop", "translation", "notes", "warnings", "refHits", "edited", "time_stamp"
)

#' TRUE when a lifted gene of this type is curated through Nextflow
#' @noRd
liftover_needs_nf <- function(type) length(type) == 1 && type %in% c("PCG", "rRNA")

#' Build a lifted annotation row and write the one-row curation inputs to `dir`
#'
#' `curate_params` is the parsed list or the JSON text stored in curate_opts.
#' Returns the state list used by the curation runners and liftover_finish().
#' @noRd
liftover_prepare <- function(hit, ref_feature, accession, unit, assembly, coverage_fn,
                             curate_params, genetic_code, ref_dir, template,
                             max_blast_hits = 10, blast_ref_files = character(),
                             dir = tempfile("liftover_")) {
  json <- if (is.character(curate_params)) curate_params else jsonlite::toJSON(curate_params)
  if (is.character(curate_params)) curate_params <- jsonlite::fromJSON(curate_params)
  asm_s <- assembly[[1]]
  L <- length(asm_s)
  contig <- sub("\\s.*", "", names(assembly)[1])
  type <- ref_feature$type[1]
  gene <- ref_feature$gene[1]
  row <- data.frame(
    contig = contig, type = type, gene = gene,
    product = liftover_product(gene, type, template),
    pos1 = hit$pos1[1], pos2 = hit$pos2[1],
    length = as.integer(circ_len(hit$pos1[1], hit$pos2[1], L)),
    direction = hit$direction[1], anticodon = NA_character_,
    tool = paste0(accession, "_liftover"),
    start_codon = NA_character_, stop_codon = NA_character_,
    partial_start = 0L, partial_stop = 0L, translation = NA_character_,
    notes = sprintf("Liftover from %s: %.1f%% identity, %.1f%% gaps", accession, hit$pctid[1], hit$gap_pct[1]),
    warnings = NA_character_, refHits = NA_character_
  )
  if (type == "PCG") {
    cds <- as.character(feature_nt(asm_s, row$pos1, row$pos2, row$direction))
    n <- nchar(cds)
    k <- if (n %% 3L == 0L) 3L else n %% 3L
    row$start_codon <- substr(cds, 1L, 3L)
    row$stop_codon <- substr(cds, n - k + 1L, n)
    row$translation <- as.character(Biostrings::translate(
      Biostrings::DNAString(substr(cds, 1L, n - k)),
      genetic.code = Biostrings::getGeneticCode(as.character(genetic_code)),
      if.fuzzy.codon = "solve"
    ))
  }
  # Inert row at the contig end: keeps one-row codon columns from being read as
  # logical (a lone "T") and gives validation the contig length when there is no coverage.
  pad <- row
  pad[1, ] <- NA
  pad[1, c("contig", "type", "gene", "pos1", "pos2", "length", "start_codon", "stop_codon", "translation")] <-
    list(contig, "ctrl", ".liftover_pad", L, L, 1L, "ATG", "TAA", "X")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir <- normalizePath(dir)
  ann_fn <- file.path(dir, "liftover_annotations_1.csv")
  asm_fn <- file.path(dir, "liftover_assembly_1.fasta")
  utils::write.csv(rbind(row, pad), ann_fn, row.names = FALSE)
  Biostrings::writeXStringSet(assembly, asm_fn)
  list(
    row = row, pad = pad, type = type, unit = unit, assembly = assembly, dir = dir,
    ann_fn = ann_fn, asm_fn = asm_fn, coverage_fn = coverage_fn, ref_dir = ref_dir,
    curate_params = curate_params, params_b64 = base64enc::base64encode(charToRaw(as.character(json))),
    genetic_code = genetic_code, max_blast_hits = max_blast_hits, blast_ref_files = blast_ref_files
  )
}

#' Curate the prepared row in this R process
#'
#' Writes the curated table and assembly under `<dir>/curate`; returns NULL on
#' success, else the reason curation was skipped.
#' @noRd
liftover_curate_local <- function(state) {
  st <- state
  if (st$type == "PCG" && !isTRUE(dir.exists(file.path(st$ref_dir, "featureProt")))) {
    return("curation reference database not found")
  }
  out_dir <- file.path(st$dir, "curate")
  dir.create(out_dir, showWarnings = FALSE)
  reason <- NULL
  cur <- tryCatch(
    withCallingHandlers(
      curate_mito_core(
        annotations_fn = st$ann_fn, assembly_fn = st$asm_fn, coverage_fn = st$coverage_fn,
        genetic_code = st$genetic_code, out_dir = out_dir, max_blast_hits = st$max_blast_hits,
        params = st$curate_params, ref_dir = st$ref_dir, feature_trim = FALSE,
        flip_rRNA_minus_strand = FALSE, ref_based_rc = FALSE
      ),
      warning = function(w) {
        if (grepl("had status", conditionMessage(w))) stop(conditionMessage(w), call. = FALSE)
      }
    ),
    error = function(e) {
      while (inherits(e$parent, "condition")) e <- e$parent
      reason <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(cur)) return(reason)
  utils::write.csv(as.data.frame(cur), file.path(out_dir, basename(st$ann_fn)), row.names = FALSE)
  NULL
}

#' Finish a lift: guard, validate, strip the pad row
#'
#' `curated_dir` holds `curate/` (and optionally `validate/`) outputs for the
#' prepared files; NULL adds the row uncurated with `skip_reason`. Returns
#' list(row = 1-row data.frame with the annotations table columns, curated).
#' Stops with a user-readable message on failure.
#' @noRd
liftover_finish <- function(state, curated_dir, skip_reason = NULL) {
  st <- state
  curated <- !is.null(curated_dir)
  if (curated) {
    out_asm <- Biostrings::readDNAStringSet(file.path(curated_dir, "curate", basename(st$asm_fn)))
    if (!identical(as.character(out_asm), as.character(st$assembly))) {
      stop("Liftover aborted: curation changed the assembly", call. = FALSE)
    }
    cur <- utils::read.csv(file.path(curated_dir, "curate", basename(st$ann_fn)), na.strings = c("", "NA"))
  } else {
    cur <- rbind(st$row, st$pad)
    reason <- trimws(gsub("[;[:space:]]+", " ", skip_reason %||% ""))
    if (!nzchar(reason)) reason <- "unknown error"
    cur$warnings[1] <- paste0("Liftover: curation skipped (", reason, ")")
  }
  keep <- which(cur$gene != ".liftover_pad")
  if (length(keep) != 1L || is.na(cur$pos1[keep]) || is.na(cur$pos2[keep])) {
    stop("Liftover failed: curation removed the annotation.", call. = FALSE)
  }
  tsv <- file.path(curated_dir %||% "", "validate", sub("\\.csv$", ".tsv", basename(st$ann_fn)))
  if (curated && file.exists(tsv)) {
    val <- utils::read.delim(tsv, quote = "", na.strings = "", colClasses = "character")
  } else {
    vdir <- file.path(st$dir, "validate")
    dir.create(vdir, showWarnings = FALSE)
    val_fn <- file.path(vdir, basename(st$ann_fn))
    utils::write.csv(cur, val_fn, row.names = FALSE)
    val <- validate_mito_core(
      annotations_fn = val_fn, coverage_fn = st$coverage_fn,
      params = st$curate_params, out_dir = vdir
    )$annotations
  }
  out <- cur[keep, , drop = FALSE]
  out$warnings <- val$warnings[val$gene != ".liftover_pad"]
  if (!curated) {
    # no reference hits were fetched, so these checks cannot be judged
    w <- trimws(strsplit(out$warnings %|NA|% "", ";", fixed = TRUE)[[1]])
    w <- w[!grepl("^(low reference similarity|check reference (start|stop) alignment)$", w)]
    out$warnings <- paste(w, collapse = "; ")
  }
  out$ID <- st$unit$ID
  out$path <- st$unit$path
  out$scaffold <- st$unit$scaffold
  out$edited <- 1
  out$time_stamp <- as.numeric(Sys.time())
  out <- out[, LIFTOVER_ANN_COLS]
  rownames(out) <- NULL
  list(row = out, curated = curated)
}

#' Build a lifted annotation row, curate and validate it alone in this process
#'
#' Returns list(row = 1-row data.frame with the annotations table columns,
#' curated = TRUE/FALSE). Stops with a user-readable message on failure.
#' @noRd
liftover_annotate <- function(hit, ref_feature, accession, unit, assembly, coverage_fn,
                              curate_params, genetic_code, ref_dir, template) {
  st <- liftover_prepare(hit, ref_feature, accession, unit, assembly, coverage_fn,
                         curate_params, genetic_code, ref_dir, template)
  on.exit(unlink(st$dir, recursive = TRUE), add = TRUE)
  reason <- liftover_curate_local(st)
  liftover_finish(st, if (is.null(reason)) st$dir, reason)
}

#' Candidate-reference BLAST JSONs a pipeline run gives curation for a sample
#' @noRd
liftover_blast_ref_files <- function(dir_out, id, assemble_opts) {
  d <- file.path(dir_out, id, "assemble", assemble_opts)
  f <- sort(Sys.glob(file.path(d, "blast_ref_*", "remote_blast_ref.json")))
  if (length(f) == 0 && file.exists(file.path(d, "remote_blast_ref.json"))) f <- file.path(d, "remote_blast_ref.json")
  f
}

#' Download (when `url` is set) and unpack a curation DB tarball into `out`
#'
#' Base R only, so the app can run it in a child R process.
#' @noRd
liftover_db_fetch <- function(url, tb, staging, out, clade) {
  op <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(op), add = TRUE)
  if (nzchar(url)) {
    ok <- tryCatch(
      utils::download.file(url, tb, quiet = TRUE, mode = "wb") == 0,
      error = function(e) FALSE, warning = function(w) FALSE
    )
    if (!ok) stop("Could not download the curation reference database (", url, ").", call. = FALSE)
  }
  if (utils::untar(tb, exdir = staging) != 0) {
    stop("Could not unpack the curation reference database (", tb, ").", call. = FALSE)
  }
  src <- file.path(staging, clade)
  if (!dir.exists(file.path(src, "featureProt"))) src <- staging
  if (!dir.exists(file.path(src, "featureProt"))) {
    stop("The curation reference database (", basename(tb), ") has no featureProt folder.", call. = FALSE)
  }
  if (!suppressWarnings(file.rename(src, out)) && !dir.exists(file.path(out, "featureProt"))) {
    stop("Could not cache the curation reference database in ", out, call. = FALSE)
  }
  invisible(out)
}

#' Local directory holding the curation reference DB for a unit
#'
#' A local `<ref_dir>/<ref_db>` directory is used as is. A local tarball or a
#' URL is unpacked once into `<project>/.liftover/ref_dbs/<ref_db>` and reused.
#' With `background = TRUE` a missing cache is fetched by a child R process,
#' which is returned instead of the path; call again once it exits with 0.
#' @noRd
liftover_ref_db_dir <- function(ref_dir, ref_db, project_dir, background = FALSE) {
  clade <- sub("\\.tar\\.gz$", "", basename(ref_db %|NA|% ""))
  if (!nzchar(ref_dir %|NA|% "") || !nzchar(clade)) {
    stop("No curation reference database is set for this assembly.", call. = FALSE)
  }
  remote <- grepl("^(https?|ftp)://", ref_dir)
  if (!remote) {
    d <- file.path(ref_dir, clade)
    if (dir.exists(d)) return(normalizePath(d))
  }
  out <- file.path(project_dir, ".liftover", "ref_dbs", clade)
  if (dir.exists(out)) return(normalizePath(out))
  tb <- paste0(file.path(ref_dir, clade), ".tar.gz")
  if (!remote && !file.exists(tb)) {
    stop("Curation reference database not found: ", file.path(ref_dir, clade), call. = FALSE)
  }
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(".unpack_", tmpdir = dirname(out))
  dir.create(file.path(tmp, "x"), recursive = TRUE)
  url <- if (remote) tb else ""
  if (remote) tb <- file.path(tmp, basename(url))
  if (!background) {
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    liftover_db_fetch(url, tb, file.path(tmp, "x"), out, clade)
    return(normalizePath(out))
  }
  fn <- paste(deparse(liftover_db_fetch), collapse = "\n")
  call <- sprintf(
    "r <- tryCatch(f(%s, %s, %s, %s, %s), error = function(e) e)",
    deparse(url), deparse(tb), deparse(file.path(tmp, "x")), deparse(out), deparse(clade)
  )
  code <- paste(
    paste("f <-", fn), call, sprintf("unlink(%s, recursive = TRUE)", deparse(tmp)),
    "if (inherits(r, 'error')) { cat(conditionMessage(r)); quit(status = 1) }",
    sep = "\n"
  )
  p <- processx::process$new(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", "-e", code),
    stdout = "|", stderr = "2>&1", cleanup_tree = TRUE
  )
  attr(p, "staging") <- tmp
  p
}

#' Nextflow config added after the project's: liftover runs on this machine
#' @noRd
LIFTOVER_NF_CONFIG <- "process {
  withName: 'liftover_curate' {
    executor = 'local'
    cpus = 1
    memory = '2 GB'
  }
}"

#' Command, arguments, working directory, and environment for a liftover run
#'
#' Launcher, NXF_VER pin, and native environment follow the app's pipeline run.
#' @noRd
liftover_nf_command <- function(state, db_dir, project_dir, config = file.path(project_dir, ".config")) {
  st <- state
  d <- st$dir
  nat <- read_config_executor(config)$native_activate
  nat_env <- if (!is.null(nat)) native_env(nat) else character()
  pin <- if (length(nat_env)) native_nf_pin(nat_env) else nf_pin_version()
  cov <- if (is.null(st$coverage_fn)) app_sys("nextflow", "assets", "NO_FILE") else normalizePath(st$coverage_fn)
  cache <- liftover_sif_cache_env(config, project_dir, nat_env)
  refs <- if (length(st$blast_ref_files)) normalizePath(st$blast_ref_files) else app_sys("nextflow", "modules", "empty_remote_blast_ref.json")
  list(
    command = "nextflow",
    args = c(
      "-log", file.path(d, "nextflow.log"), "run", app_sys("nextflow", "liftover_curate.nf"),
      "-c", config, "-c", file.path(d, "liftover.config"), "-w", file.path(d, "work"),
      "-ansi-log", "false",
      "--lift_annotations", st$ann_fn, "--lift_assembly", st$asm_fn, "--lift_coverage", cov,
      "--lift_ref_db", db_dir, "--lift_blast_ref", paste(refs, collapse = ","), "--lift_params", st$params_b64,
      "--lift_genetic_code", as.character(st$genetic_code),
      "--lift_max_blast_hits", as.character(st$max_blast_hits), "--lift_outdir", d
    ),
    wd = d,
    env = c("current", nat_env, if (!is.na(pin)) c(NXF_VER = pin), cache)
  )
}

#' Singularity/Apptainer image cache shared with the project's pipeline runs
#'
#' `-w <scratch>/work` would move Nextflow's default cache (`<workDir>/singularity`)
#' into the per-lift folder, so point it at the pipeline's default, `<project>/work/singularity`,
#' unless the config or the environment already names a cache.
#' @noRd
liftover_sif_cache_env <- function(config, project_dir, nat_env = character()) {
  keys <- c("NXF_SINGULARITY_CACHEDIR", "NXF_APPTAINER_CACHEDIR")
  lines <- if (file.exists(config)) sub("//.*$", "", readLines(config, warn = FALSE)) else character()
  if (any(grepl("cacheDir\\s*=", lines))) return(character())
  if (any(nzchar(Sys.getenv(keys))) || any(keys %in% names(nat_env))) return(character())
  stats::setNames(rep(file.path(project_dir, "work", "singularity"), 2), keys)
}

#' Start the liftover workflow in the background; returns the processx process
#' @noRd
liftover_nf_start <- function(state, db_dir, project_dir, config = file.path(project_dir, ".config")) {
  writeLines(LIFTOVER_NF_CONFIG, file.path(state$dir, "liftover.config"))
  cmd <- liftover_nf_command(state, db_dir, project_dir, config)
  processx::process$new(
    cmd$command, cmd$args, wd = cmd$wd, env = cmd$env,
    stdout = file.path(state$dir, "nextflow.out"), stderr = "2>&1", cleanup_tree = TRUE
  )
}

#' Plain one-line reason from a failed liftover run's console output
#' @noRd
liftover_nf_reason <- function(lines) {
  lines <- trimws(lines)
  i <- match("Command error:", lines)
  if (!is.na(i)) {
    blk <- lines[-seq_len(i)]
    end <- match("", blk)
    if (!is.na(end)) blk <- blk[seq_len(end - 1L)]
    blk <- blk[nzchar(blk) & !grepl("^(Execution halted|Calls:)", blk)]
    if (length(blk)) return(blk[length(blk)])
  }
  e <- grep("^ERROR ~ ", lines, value = TRUE)
  if (length(e)) return(sub("^ERROR ~ ", "", e[1]))
  "Nextflow stopped without a message"
}

#' Display table for liftover hits: text cells, substantial-overlap genes, pass flag
#'
#' Identity is floored and gaps ceiled to one decimal so a shown value never
#' looks like it passes when it does not.
#' @noRd
liftover_hit_table <- function(hits, annotations, asm_len, min_pctid, max_gap_pct) {
  thr <- liftover_clamp(min_pctid, max_gap_pct)
  ov <- lapply(seq_len(nrow(hits)), function(i) {
    liftover_overlaps(hits$pos1[i], hits$pos2[i], annotations, asm_len)
  })
  fmt <- function(x) formatC(as.integer(x), big.mark = ",", format = "d")
  data.frame(
    position = paste0(fmt(hits$pos1), "-", fmt(hits$pos2), recycle0 = TRUE),
    strand = hits$direction,
    identity = sprintf("%.1f%%", floor(hits$pctid * 10) / 10),
    gaps = sprintf("%.1f%%", ceiling(hits$gap_pct * 10) / 10),
    overlaps = vapply(ov, function(o) {
      if (nrow(o) == 0) "none" else paste0(o$gene, " (", o$bp, " bp)", collapse = ", ")
    }, character(1)),
    clash = vapply(ov, function(o) paste(unique(o$gene[o$substantial]), collapse = ","), character(1)),
    pass = hits$pctid >= thr$min_pctid & hits$gap_pct <= thr$max_gap_pct
  )
}

#' Confirmation text for a lifted gene that clashes, or NULL when none
#' @noRd
liftover_confirm_text <- function(gene, overlap_genes, same_name) {
  and <- function(x) {
    n <- length(x)
    if (n <= 2) paste(x, collapse = " and ") else paste0(paste(x[-n], collapse = ", "), ", and ", x[n])
  }
  msg <- c(
    if (length(overlap_genes)) paste0("This placement overlaps ", and(overlap_genes), " by more than a quarter of either gene's length."),
    if (same_name) paste0("This assembly already has an annotation named ", gene, ".")
  )
  if (length(msg) == 0) return(NULL)
  paste(c(msg, paste("Add", gene, "anyway?")), collapse = " ")
}

#' Two-row alignment of a lifted placement for the dialog
#'
#' Nucleotides are the mapping alignment itself. A PCG whose placement starts
#' with a valid start codon is shown as translated protein instead.
#' @noRd
liftover_alignment <- function(hit, ref_feature, ref_seq, assembly_seq, type,
                               start_codons, sample_code, ref_code, labels) {
  nt <- list(
    aln = Biostrings::DNAStringSet(stats::setNames(c(hit$aln_gene[1], hit$aln_sample[1]), labels)),
    kind = "nt", notes = character()
  )
  if (!identical(type, "PCG")) return(nt)
  asm <- Biostrings::DNAString(liftover_clean(assembly_seq))
  cds <- as.character(feature_nt(asm, hit$pos1[1], hit$pos2[1], hit$direction[1]))
  ref <- liftover_ref_gene_seq(ref_seq, ref_feature$pos1[1], ref_feature$pos2[1], ref_feature$direction[1])
  tr <- function(x, code) {
    x <- substr(x, 1L, nchar(x) - nchar(x) %% 3L)
    aa <- Biostrings::translate(
      Biostrings::DNAString(x),
      genetic.code = Biostrings::getGeneticCode(as.character(code)), if.fuzzy.codon = "solve"
    )
    sub("\\*$", "", as.character(aa))
  }
  aa <- tr(cds, sample_code)
  n_stop <- nchar(gsub("[^*]", "", aa))
  nt$notes <- if (n_stop > 0) {
    sprintf("%d internal stop codon%s in the mapped region.", n_stop, if (n_stop == 1) "" else "s")
  } else character()
  if (length(start_codons) == 0) {
    nt$notes <- c("Start codons for this gene are not known, showing nucleotide alignment.", nt$notes)
    return(nt)
  }
  codon <- substr(cds, 1L, 3L)
  if (!codon %in% start_codons) {
    nt$notes <- c(sprintf("Placement does not start with a valid start codon (%s), showing nucleotide alignment.", codon), nt$notes)
    return(nt)
  }
  a <- pwalign::pairwiseAlignment(
    Biostrings::AAString(tr(ref, ref_code)), Biostrings::AAString(aa),
    substitutionMatrix = "BLOSUM62", type = "global"
  )
  aln <- Biostrings::AAStringSet(c(as.character(pwalign::alignedPattern(a)), as.character(pwalign::alignedSubject(a))))
  names(aln) <- labels
  list(aln = aln, kind = "aa", notes = nt$notes)
}

#' TRUE for annotations added by liftover (tool "<accession>_liftover")
#' @noRd
is_liftover_row <- function(tool) isTRUE(endsWith(as.character(tool), "_liftover"))

# A 5'-partial CDS whose frame starts at base 2 or 3 must keep that offset in
# the database, be left alone by curation, and export with codon_start.

cs_db <- function() {
  td <- withr::local_tempdir(.local_envir = parent.frame())
  mapping <- file.path(td, "mapping.csv")
  utils::write.csv(
    data.frame(ID = "S1", Taxon = "Danio rerio",
               R1 = "S1_R1.fastq.gz", R2 = "S1_R2.fastq.gz"),
    mapping, row.names = FALSE
  )
  suppressMessages(new_db(db_path = file.path(td, ".sqlite"), mapping_fn = mapping))
  writeLines(c("params {", "    rawDir = 'raw'", "}", "process {",
               "    executor = 'local'", "}"), file.path(td, ".config"))
  td
}

test_that("new DBs have annotations.codon_start and old DBs gain it", {
  td <- cs_db()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  expect_true("codon_start" %in% DBI::dbListFields(con, "annotations"))
  if ("codon_start" %in% DBI::dbListFields(con, "annotations")) {
    DBI::dbExecute(con, "ALTER TABLE annotations DROP COLUMN codon_start")
  }
  DBI::dbDisconnect(con)

  suppressMessages(backwards_compatibility(path = td, update_config = FALSE))

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  expect_true("codon_start" %in% DBI::dbListFields(con, "annotations"))
})

test_that("validate.nf writes codon_start with time_stamp and edited after it", {
  nf <- readLines(system.file("nextflow", "modules", "validate.nf", package = "MitoPilot"))
  src <- paste(nf, collapse = "\n")
  grab <- function(name) {
    m <- regmatches(src, regexpr(paste0("def ", name, " = \\[[^]]*\\]"), src))
    eval(parse(text = sub("^def \\w+ = \\[", "c(", sub("\\]$", ")", m))))
  }
  ins <- grab("insCols")
  tsv <- grab("tsvCols")
  expect_true("codon_start" %in% ins)
  expect_true("codon_start" %in% tsv)
  # TSV columns bind at 4 + k (k from 0), right after ID/path/scaffold
  expect_identical(ins[4:(3 + length(tsv))], tsv)
  ts_idx <- as.integer(sub(".*setStr\\(ins, (\\d+), ts\\).*", "\\1",
                           grep("setStr\\(ins, \\d+, ts\\)", nf, value = TRUE)))
  ed_idx <- as.integer(sub(".*ins\\.setInt\\((\\d+), 0\\).*", "\\1",
                           grep("ins\\.setInt\\(\\d+, 0\\)", nf, value = TRUE)))
  expect_equal(ins[ts_idx], "time_stamp")
  expect_equal(ins[ed_idx], "edited")
  expect_equal(length(ins), ed_idx)
})

test_that("curate boundary loop skips codon_start > 1 rows", {
  notes <- c(NA, NA, NA, "JOIN: mode=exon group=1", NA)
  cs <- c(NA, 1L, 2L, 1L, 3L)
  expect_identical(skip_boundary_curation(notes, cs),
                   c(FALSE, FALSE, TRUE, TRUE, TRUE))
  # tool path has no codon_start column
  expect_identical(skip_boundary_curation(notes, NULL),
                   c(FALSE, FALSE, FALSE, TRUE, FALSE))
})

test_that("splice_join_cds drops the 5' member's offset before translating", {
  # ATG GCC AAA on + strand, preceded by 2 junk bases in the 5' exon
  seq <- Biostrings::DNAString("CCATGGCCAAATAA")
  members <- data.frame(
    pos1 = c(1L, 9L), pos2 = c(8L, 14L), direction = "+",
    start_codon = c("CCA", NA), stop_codon = c(NA, "TAA"),
    partial_start = c(1L, 0L), partial_stop = c(0L, 0L),
    codon_start = c(3L, 1L), stringsAsFactors = FALSE
  )
  sp <- splice_join_cds(members, seq, Biostrings::getGeneticCode("2"))
  expect_equal(sp$translation, "MAK")
  expect_equal(sp$codon_start, 3L)
  members$codon_start <- NULL
  expect_equal(splice_join_cds(members, seq, Biostrings::getGeneticCode("2"))$codon_start, 1L)
})

L <- 3000L
withr::with_seed(4, {
  cs_seq <- paste(sample(c("A", "C", "G", "T"), L, replace = TRUE), collapse = "")
})

cs_row <- function(gene, pos1, pos2, codon_start = NA_integer_, partial_start = 0L) {
  data.frame(
    ID = "s1", path = 1L, scaffold = 1L, contig = "s1.1.1", type = "PCG",
    gene = gene, product = gene, pos1 = as.integer(pos1), pos2 = as.integer(pos2),
    length = as.integer(pos2 - pos1 + 1L), direction = "+",
    start_codon = "ATG", stop_codon = "TAA",
    translation = strrep("M", 30L), anticodon = NA_character_,
    partial_start = partial_start, partial_stop = 0L, codon_start = codon_start,
    notes = NA_character_, refHits = "{}", warnings = NA_character_,
    stringsAsFactors = FALSE
  )
}

cs_export <- function(annotations) {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  out_dir <- file.path(d, "out")
  dir.create(out_dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  DBI::dbWriteTable(con, "annotations", annotations)
  DBI::dbWriteTable(con, "assemblies", data.frame(
    ID = "s1", path = 1L, scaffold = 1L, ignore = 0L, topology = "circular",
    sequence = cs_seq, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "export", data.frame(
    ID = "s1", path = 1L, scaffold = 1L, export_group = "g1",
    export_time_stamp = NA_integer_, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "annotate", data.frame(
    ID = "s1", path = 1L, scaffold = 1L, topology = "circular", partial = "no",
    curate_opts = "default", stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "curate_opts", data.frame(
    curate_opts = "default",
    params = '{"default_rules":{"PCG":{"start_codons":["ATG"],"stop_codons":["TAA"]}}}',
    linear_complete = 0L,
    stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "samples", data.frame(
    ID = "s1", Taxon = "Testus testus", genetic_code = 2L, stringsAsFactors = FALSE))
  DBI::dbWriteTable(con, "assemble", data.frame(
    ID = "s1", blast_accession = "NC_000001", blast_accession_auto = 0L,
    poor_blast_ref = "ok", stringsAsFactors = FALSE))
  DBI::dbDisconnect(con)

  suppressMessages(suppressWarnings(export_files(
    group = "g1", out_dir = out_dir, generateAAalignments = FALSE,
    gene_export = FALSE, review = FALSE, summary_csv = FALSE
  )))
  structure(readLines(file.path(out_dir, "s1", "export", "s1.tbl")),
            gff = readLines(file.path(out_dir, "s1", "export", "s1.gff")))
}

test_that("export writes codon_start 2 for an offset partial CDS", {
  tbl <- cs_export(rbind(
    cs_row("cox1", 101, 400, codon_start = 2L, partial_start = 1L),
    cs_row("cox2", 501, 800)
  ))
  expect_true(any(grepl("^<101\t400\tgene$", tbl)))
  expect_true(any(grepl("^\t\t\tcodon_start\t\\s*2$", tbl)))
  # a complete gene with a valid start gets no codon_start line
  expect_equal(sum(grepl("codon_start", tbl)), 1L)
  cds <- strsplit(grep("\tCDS\t", attr(tbl, "gff"), value = TRUE), "\t")
  phase <- setNames(vapply(cds, `[`, "", 8L), vapply(cds, `[`, "", 4L))
  expect_equal(phase[["101"]], "1")
  expect_equal(phase[["501"]], "0")
})

test_that("GFF phase follows codon_start across JOIN exons", {
  a <- rbind(cs_row("nad5", 101, 201, codon_start = 2L, partial_start = 1L),
             cs_row("nad5", 301, 400))
  a$notes <- "JOIN: mode=exon group=1"
  cds <- strsplit(grep("\tCDS\t", attr(cs_export(a), "gff"), value = TRUE), "\t")
  # 101 bp exon read from base 2: 100 bases, 1 left over, so exon 2 starts at phase 2
  expect_equal(vapply(cds, `[`, "", 8L), c("1", "2"))
})

test_that("shipped user annotation fixtures parse", {
  td <- system.file("test_data", package = "MitoPilot")
  map <- utils::read.csv(file.path(td, "mapping_test_userAnnot.csv"))
  expect_setequal(list.files(file.path(td, "annotations")), map$Annotation)
  unmapped <- list(
    GB_NCBI = character(0), GB_TBL = character(0), GB_BED = character(0),
    GB_TYPOS = c("COX11", "NDH4K", "ATP88", "cytochrome oxidase subnit I", "tRNA-Xxx",
                 "rrnX", "hypothetical protein", "ORF123"),
    GB_DROP = c("", "control region", "putative light strand origin of replication",
                "tandem repeat"),
    GB_MITOS = c("OH", "OL"))
  p <- list()
  for (i in seq_len(nrow(map))) {
    len <- unname(Biostrings::fasta.seqlengths(file.path(td, "assemblies", map$Assembly[i])))
    p[[map$ID[i]]] <- parse_user_annotation(file.path(td, "annotations", map$Annotation[i]), len)
    x <- p[[map$ID[i]]]
    g <- guess_gene(x$raw_name, x$feature_type, x$raw_product)$gene
    expect_setequal(unique(dplyr::coalesce(x$raw_name[is.na(g)], "")), unmapped[[map$ID[i]]])
  }
  ok <- c("mt-Co1" = "cox1", "MT-ND4L" = "nad4l", "cytb" = "cob", "COII" = "cox2",
          "tRNA-Leu (CUN)" = "trnL", "trnS2(tga)" = "trnS", "12S rRNA" = "rrnS",
          "l-rRNA" = "rrnL", "cytochrome c oxidase subunit III" = "cox3")
  x <- p$GB_TYPOS[match(names(ok), p$GB_TYPOS$raw_name), ]
  expect_equal(guess_gene(x$raw_name, x$feature_type, x$raw_product)$gene, unname(ok))
  x <- p$GB_MITOS[match(c("cox1", "nad4l", "trnL1(tag)", "trnS2(tga)"), p$GB_MITOS$raw_name), ]
  expect_equal(guess_gene(x$raw_name, x$feature_type, x$raw_product)$gene,
               c("cox1", "nad4l", "trnL", "trnS"))
  d <- p$GB_DROP
  expect_equal(d$codon_start[d$raw_name %in% "ND1"], 2L)
  expect_equal(d$partial_start[d$raw_name %in% "ND1"], 1L)
  dr <- attr(d, "dropped_types")
  expect_length(dr, 2L)
  expect_equal(dr[c("STS", "gap")], c(STS = 1L, gap = 1L))
})

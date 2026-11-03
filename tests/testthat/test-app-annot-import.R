ua_app_project <- function(td) {
  src <- app_sys("test_data/assemblies")
  asm <- file.path(td, "asm")
  ann <- file.path(td, "ann")
  dir.create(asm)
  dir.create(ann)
  file.copy(file.path(src, "UA_LINEAR.fasta"), asm)
  writeLines(c("##gff-version 3", paste("mito_contig", "src", c(
    "CDS\t100\t300\t.\t+\t0\tID=c1;Name=COX1",
    "CDS\t1000\t1200\t.\t+\t0\tID=o1;Name=ORF1"
  ), sep = "\t")), file.path(ann, "lin.gff3"))
  m <- data.frame(ID = "UA_LINEAR", Taxon = "Danio rerio", Assembly = "UA_LINEAR.fasta",
                  Topology = "linear", Annotation = "lin.gff3")
  utils::write.csv(m, file.path(td, "mapping.csv"), row.names = FALSE)
  proj <- file.path(td, "proj")
  suppressMessages(new_project_userAnnot(
    path = proj, mapping_fn = file.path(td, "mapping.csv"), assembly_path = asm,
    annotation_path = ann, no_raw_data = TRUE, Rproj = FALSE, executor = "local"))
  proj
}

ua_app_project2 <- function(td) {
  src <- app_sys("test_data/assemblies")
  asm <- file.path(td, "asm")
  ann <- file.path(td, "ann")
  dir.create(asm)
  dir.create(ann)
  file.copy(file.path(src, c("UA_LINEAR.fasta", "UA_CIRCULAR.fasta")), asm)
  writeLines(c("##gff-version 3", paste("mito_contig", "src", c(
    "CDS\t100\t300\t.\t+\t0\tID=c1;Name=COX1",
    "CDS\t500\t700\t.\t-\t0\tID=c2;Name=COX1",
    "CDS\t1000\t1200\t.\t+\t0\tID=o1;Name=ORF1"
  ), sep = "\t")), file.path(ann, "lin.gff3"))
  writeLines(c("mito_contig\t49\t250\tCOX1\t0\t+", "mito_contig\t899\t1100\tmystery\t0\t+"),
             file.path(ann, "circ.bed"))
  m <- data.frame(ID = c("UA_LINEAR", "UA_CIRCULAR"), Taxon = "Danio rerio",
                  Assembly = c("UA_LINEAR.fasta", "UA_CIRCULAR.fasta"),
                  Topology = c("linear", "circular"), Annotation = c("lin.gff3", "circ.bed"))
  utils::write.csv(m, file.path(td, "mapping.csv"), row.names = FALSE)
  proj <- file.path(td, "proj")
  suppressMessages(new_project_userAnnot(
    path = proj, mapping_fn = file.path(td, "mapping.csv"), assembly_path = asm,
    annotation_path = ann, no_raw_data = TRUE, Rproj = FALSE, executor = "local"))
  proj
}

ua_app_session <- function(con, td) {
  ms <- shiny::MockShinySession$new()
  ms$userData$con <- con
  ms$userData$dir <- td
  ms$userData$user_annot <- TRUE
  for (f in c("refresh_assemble", "refresh_annotate", "annot_import_state")) gargoyle::init(f, session = ms)
  ms
}

test_that("run gate blocks only unconfirmed user annotation projects", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  expect_false(annot_import_blocks_run(list(userData = list(user_annot = FALSE, con = con))))
  expect_true(annot_import_blocks_run(list(userData = list(user_annot = TRUE, con = con))))
  confirm_user_annotations(con, td)
  expect_false(annot_import_blocks_run(list(userData = list(user_annot = TRUE, con = con))))
})

test_that("annot_import_notice reports unmapped counts until confirmed", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  expect_equal(annot_import_notice(con),
               "Gene names not confirmed: 1 unmapped of 2. Pipeline runs are blocked.")
  confirm_user_annotations(con, td)
  expect_null(annot_import_notice(con))
})

test_that("mp_app_header renders the notice with a Review button", {
  h <- as.character(mp_app_header(tempdir(), notice = "Gene names not confirmed"))
  expect_match(h, "mp-header-notice", fixed = TRUE)
  expect_match(h, "gene_names_header", fixed = TRUE)
  expect_false(grepl("gene_names_header", as.character(mp_app_header(tempdir())), fixed = TRUE))
})

test_that("annot_import pick, reset, and confirm", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 0")
  shiny::testServer(annot_import_server, args = list(id = "x"), session = ua_app_session(con, td), {
    rows <- function() jsonlite::fromJSON(output$tbl)$x$tag$attribs$data
    expect_equal(rows()$raw_name, c("ORF1", "COX1"))
    expect_equal(rows()$status, c("unmapped", "auto"))
    expect_equal(rows()$ftype, c("GFF3", "GFF3"))
    expect_equal(rows()$class, c("none", "PCG"))
    expect_match(output$tbl, "x-pick", fixed = TRUE)
    expect_match(output$tbl, "ctrl (control region)", fixed = TRUE)
    expect_match(as.character(output$chips$html), "Unmapped (1)", fixed = TRUE)

    st <- session$userData$annot_import_state()
    session$setInputs(pick = list(raw_name = "ORF1", feature_type = "CDS", gene = "nad6"))
    expect_gt(session$userData$annot_import_state(), st)
    gm <- DBI::dbGetQuery(con, "SELECT gene, manual FROM gene_name_map WHERE raw_name = 'ORF1'")
    expect_equal(gm$gene, "nad6")
    expect_equal(gm$manual, 1L)
    d <- annot_import_data(con, c("ORF1\rCDS", "COX1\rCDS"))
    expect_equal(d$raw_name, c("ORF1", "COX1"))
    expect_equal(d$status, c("edited", "auto"))
    expect_equal(d$class, c("PCG", "PCG"))
    expect_match(as.character(output$chips$html), "Edited (1)", fixed = TRUE)

    st <- session$userData$annot_import_state()
    session$setInputs(reset = 1)
    expect_gt(session$userData$annot_import_state(), st)
    expect_true(is.na(DBI::dbGetQuery(con, "SELECT gene FROM gene_name_map WHERE raw_name = 'ORF1'")$gene))

    session$setInputs(pick = list(raw_name = "ORF1", feature_type = "CDS", gene = "nad6"))
    st <- session$userData$annot_import_state()
    session$setInputs(confirm = 1)
    expect_gt(session$userData$annot_import_state(), st)
    expect_true(user_annot_confirmed(con))
    expect_true(file.exists(file.path(td, "user_annotations", "UA_LINEAR_annotations.csv")))
    expect_equal(DBI::dbGetQuery(con, "SELECT annotate_switch FROM annotate")$annotate_switch, 1L)
  })
})

test_that("gene_name_map_samples lists each sample's features for a name", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project2(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  s <- gene_name_map_samples(con, "COX1", "CDS")
  expect_equal(s$ID, c("UA_LINEAR", "UA_LINEAR"))
  expect_equal(s$pos1, c(100L, 500L))
  expect_equal(s$n_features, c(2L, 2L))
  expect_equal(s$direction, c("+", "-"))
  expect_equal(s$source_file, c("lin.gff3", "lin.gff3"))
  expect_true(all(c("feature_type", "partial_start", "partial_stop", "codon_start",
                    "segment_group", "pseudo", "note") %in% names(s)))
  b <- gene_name_map_samples(con, "COX1", "")
  expect_equal(b$ID, "UA_CIRCULAR")
  expect_equal(b$pos1, 50L)
  expect_equal(nrow(gene_name_map_samples(con, "nope", "CDS")), 0L)
})

test_that("annot_import_counts matches the gene_name_map and the table", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project2(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  n <- annot_import_counts(con)
  gm <- DBI::dbGetQuery(con, "SELECT gene, manual FROM gene_name_map")
  expect_equal(n$all, nrow(gm))
  expect_equal(n$unmapped, sum(is.na(gm$gene)))
  expect_equal(n$edited, 0L)
  expect_equal(n$auto, sum(!is.na(gm$gene)))
  set_gene_name_map(con, "mystery", "", "nad1")
  n <- annot_import_counts(con)
  d <- annot_import_data(con)
  expect_equal(n$edited, 1L)
  expect_equal(c(n$unmapped, n$edited, n$auto),
               as.integer(table(factor(d$status, c("unmapped", "edited", "auto")))))
  expect_equal(d$ftype[d$raw_name == "COX1" & d$feature_type == ""], "BED")
  expect_match(d$ids[d$raw_name == "COX1" & d$feature_type == ""], "UA_CIRCULAR", fixed = TRUE)
})

test_that("annot_import pick works for BED names", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project2(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  shiny::testServer(annot_import_server, args = list(id = "x"), session = ua_app_session(con, td), {
    expect_true("mystery" %in% jsonlite::fromJSON(output$tbl)$x$tag$attribs$data$raw_name)
    session$setInputs(pick = list(raw_name = "mystery", feature_type = "", gene = "trnF"))
    gm <- DBI::dbGetQuery(con, "SELECT gene, manual FROM gene_name_map WHERE raw_name = 'mystery'")
    expect_equal(gm$gene, "trnF")
    expect_equal(gm$manual, 1L)
    d <- annot_import_data(con)
    expect_equal(d$class[d$raw_name == "mystery"], "tRNA")
  })
})

test_that("unnamed rows get a label and the table builds without detail queries", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project2(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE user_annotations SET raw_name = '' WHERE raw_name = 'mystery'")
  DBI::dbExecute(con, "UPDATE gene_name_map SET raw_name = '' WHERE raw_name = 'mystery'")
  DBI::dbExecute(con, "INSERT INTO gene_name_map (raw_name, feature_type, n_samples, n_features)
                       VALUES ('', 'D-loop', 1, 1)")
  testthat::local_mocked_bindings(gene_name_map_samples = function(...) stop("eager detail query"))
  html <- as.character(htmltools::as.tags(annot_import_table(con, shiny::NS("x"))))
  expect_match(html, "(unnamed BED feature)", fixed = TRUE)
  expect_match(html, "(unnamed D-loop)", fixed = TRUE)
  expect_match(html, "mpGN.detail", fixed = TRUE)
  d <- annot_import_data(con)
  expect_equal(d$raw_name[d$label == "(unnamed BED feature)"], "")
})

test_that("annot_import open sends a detail table and adds the name index", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_app_project2(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DROP INDEX IF EXISTS idx_user_annotations_name")
  idx <- function() DBI::dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type = 'index'
                                          AND name = 'idx_user_annotations_name'")$name
  expect_length(idx(), 0)
  sent <- NULL
  ms <- ua_app_session(con, td)
  ms$sendCustomMessage <- function(type, message) if (type == "mp-gn-detail") sent <<- message
  shiny::testServer(annot_import_server, args = list(id = "x"), session = ms, {
    session$setInputs(open = list(raw_name = "COX1", feature_type = "CDS", id = "x-tbl-d-3"))
    expect_equal(sent$id, "x-tbl-d-3")
    expect_match(sent$html, "<table", fixed = TRUE)
    expect_equal(lengths(regmatches(sent$html, gregexpr("UA_LINEAR", sent$html))), 2L)
    expect_equal(idx(), "idx_user_annotations_name")
  })
  expect_match(as.character(annot_import_detail(gene_name_map_samples(con, "COX1", ""))), "UA_CIRCULAR",
               fixed = TRUE)
})

test_that(".ua_note_name cuts long notes at a word boundary", {
  expect_equal(.ua_note_name("  short note "), "short note")
  long <- "putative control region with tandem repeats and a long tail"
  out <- .ua_note_name(long)
  expect_equal(out, "putative control region with tandem...")
  expect_lte(nchar(out), 43L)
  expect_equal(.ua_note_name(strrep("x", 50)), paste0(strrep("x", 40), "..."))
  expect_true(is.na(.ua_note_name(NA_character_)))
})

test_that("export metadata picker skips the Annotation mapping column", {
  expect_equal(export_metadata_cols(c("ID", "Annotation", "sex"), character(0)), "sex")
})

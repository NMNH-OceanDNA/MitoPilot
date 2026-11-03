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

ua_app_session <- function(con, td) {
  ms <- shiny::MockShinySession$new()
  ms$userData$con <- con
  ms$userData$dir <- td
  ms$userData$user_annot <- TRUE
  for (f in c("refresh_assemble", "refresh_annotate")) gargoyle::init(f, session = ms)
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
    expect_match(output$tbl, "x-pick", fixed = TRUE)

    session$setInputs(pick = list(raw_name = "ORF1", feature_type = "CDS", gene = "nad6"))
    gm <- DBI::dbGetQuery(con, "SELECT gene, manual FROM gene_name_map WHERE raw_name = 'ORF1'")
    expect_equal(gm$gene, "nad6")
    expect_equal(gm$manual, 1L)
    expect_equal(rows()$status[rows()$raw_name == "ORF1"], "edited")

    session$setInputs(reset = 1)
    expect_true(is.na(DBI::dbGetQuery(con, "SELECT gene FROM gene_name_map WHERE raw_name = 'ORF1'")$gene))

    session$setInputs(pick = list(raw_name = "ORF1", feature_type = "CDS", gene = "nad6"))
    session$setInputs(confirm = 1)
    expect_true(user_annot_confirmed(con))
    expect_true(file.exists(file.path(td, "user_annotations", "UA_LINEAR_annotations.csv")))
    expect_equal(DBI::dbGetQuery(con, "SELECT annotate_switch FROM annotate")$annotate_switch, 1L)
  })
})

test_that("export metadata picker skips the Annotation mapping column", {
  expect_equal(export_metadata_cols(c("ID", "Annotation", "sex"), character(0)), "sex")
})

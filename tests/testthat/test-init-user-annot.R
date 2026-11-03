ua_fixture <- function(td, lin_gff = NULL, circ_gff = NULL, multi = FALSE,
                       lin_ext = "gff3") {
  src <- app_sys("test_data/assemblies")
  asm_dir <- file.path(td, "asm")
  ann_dir <- file.path(td, "ann")
  dir.create(asm_dir)
  dir.create(ann_dir)
  lin_fa <- if (multi) "UA_MULTI_ONE.fasta" else "UA_LINEAR.fasta"
  file.copy(file.path(src, c(lin_fa, "UA_CIRCULAR.fasta")), asm_dir)
  sid <- sub("\\s.*$", "", names(Biostrings::readDNAStringSet(file.path(asm_dir, lin_fa))))[1]
  gff <- function(rows) c("##gff-version 3", paste(sid, "src", rows, sep = "\t"))
  lin_gff <- lin_gff %||% c(
    "CDS\t100\t300\t.\t+\t0\tID=c1;Name=COI",
    "CDS\t400\t600\t.\t+\t0\tID=c2;Name=COI CDS 5' Partial CDS",
    "tRNA\t700\t770\t.\t-\t.\tID=t1;Name=tRNA-Phe(gaa)",
    "repeat_region\t800\t850\t.\t+\t.\tID=r1",
    "microsatellite\t900\t950\t.\t+\t.\tID=m1"
  )
  writeLines(gff(lin_gff), file.path(ann_dir, paste0("lin.", lin_ext)))
  circ_gff <- circ_gff %||% c("region\t1\t17265\t.\t+\t.\tID=reg;Is_circular=true",
                              "CDS\t17000\t17400\t.\t+\t0\tID=c1;Name=COI")
  writeLines(gff(circ_gff), file.path(ann_dir, "circ.gff"))
  m <- data.frame(ID = c("UA_LINEAR", "UA_CIRCULAR"), Taxon = "Danio rerio",
                  Assembly = c(lin_fa, "UA_CIRCULAR.fasta"),
                  Topology = c("linear", "circular"),
                  Annotation = c(paste0("lin.", toupper(lin_ext)), "circ.gff"))
  file.rename(file.path(ann_dir, paste0("lin.", lin_ext)),
              file.path(ann_dir, paste0("lin.", toupper(lin_ext))))
  map_csv <- file.path(td, "mapping.csv")
  utils::write.csv(m, map_csv, row.names = FALSE)
  list(map = map_csv, asm = asm_dir, ann = ann_dir, proj = file.path(td, "proj"))
}

ua_init <- function(fx, ...) {
  suppressMessages(new_project_userAnnot(
    path = fx$proj, mapping_fn = fx$map, assembly_path = fx$asm,
    annotation_path = fx$ann, no_raw_data = TRUE, Rproj = FALSE, executor = "local", ...))
}

test_that("new_project_userAnnot builds tables, guesses, and seeds opts", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  fx <- ua_fixture(withr::local_tempdir())
  msgs <- capture_messages(new_project_userAnnot(
    path = fx$proj, mapping_fn = fx$map, assembly_path = fx$asm,
    annotation_path = fx$ann, no_raw_data = TRUE, Rproj = FALSE, executor = "local"))
  expect_true(any(grepl("review gene names", msgs)))
  expect_true(any(grepl("(unmapped)", msgs, fixed = TRUE)))
  expect_true(any(grepl("UA_LINEAR: 4 kept, 1 dropped (microsatellite 1)", msgs, fixed = TRUE)))
  td <- fx$proj
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  expect_true(is_user_annot(con))
  gm <- DBI::dbReadTable(con, "gene_name_map")
  expect_equal(gm$gene[gm$raw_name == "COI"], "cox1")
  expect_equal(gm$n_samples[gm$raw_name == "COI"], 2L)
  expect_true(all(gm$confirmed == 0L))
  expect_true(all(gm$manual == 0L))
  ua <- DBI::dbReadTable(con, "user_annotations")
  lin <- ua[ua$ID == "UA_LINEAR", ]
  expect_equal(lin$contig[1], "UA_LINEAR.1.1")
  expect_equal(lin$partial_start[lin$raw_name == "COI CDS 5' Partial CDS"], 1L)
  expect_equal(lin$anticodon[lin$raw_name == "tRNA-Phe(gaa)"], "GAA")
  expect_equal(nrow(lin), 4L)
  circ <- ua[ua$ID == "UA_CIRCULAR", ]
  expect_true(circ$pos1 > circ$pos2)
  expect_equal(circ$contig_length, 17265L)
  expect_true(file.exists(file.path(td, "user_annotations", "raw", "UA_LINEAR.gff3")))
  expect_true(file.exists(file.path(td, "user_annotations", "raw", "UA_CIRCULAR.gff")))
  ao <- DBI::dbReadTable(con, "annotate_opts")
  expect_equal(ao$use_mitos, 0L)
  expect_equal(ao$use_trnaScan, 0L)
  expect_equal(ao$coverage_trim, 0L)
  expect_true(is.na(ao$start_gene))
  ass <- DBI::dbReadTable(con, "assemble_opts")
  expect_equal(ass$join_scaffolds, 0L)
  expect_true(all(DBI::dbReadTable(con, "find_mito_opts")$find_mitogenome == 0L))
  expect_true(file.exists(file.path(td, ".config")))
})

test_that("BED rows key on an empty feature_type", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- withr::local_tempdir()
  fx <- ua_fixture(td)
  writeLines(c("mito_contig\t99\t300\tCOI\t0\t+", "mito_contig\t500\t600\tmystery\t0\t+"),
             file.path(fx$ann, "lin.bed"))
  m <- utils::read.csv(fx$map)
  m$Annotation[1] <- "lin.bed"
  utils::write.csv(m, fx$map, row.names = FALSE)
  ua_init(fx)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(fx$proj, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  gm <- DBI::dbReadTable(con, "gene_name_map")
  expect_equal(gm$feature_type[gm$raw_name == "COI" & gm$feature_type == ""], "")
  expect_true(is.na(gm$gene[gm$raw_name == "mystery"]))
  expect_true(file.exists(file.path(fx$proj, "user_annotations", "raw", "UA_LINEAR.bed")))
})

test_that("init errors: multi-contig FASTA, seqid mismatch, wrap on linear, past end, no features", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  fx <- ua_fixture(withr::local_tempdir(), multi = TRUE)
  err <- tryCatch(ua_init(fx), error = conditionMessage)
  expect_match(err, "one sequence")
  expect_no_match(err, "find_mitogenome")
  expect_false(dir.exists(fx$proj))
  expect_length(gregexpr("holds 101 sequences", err)[[1]], 1L)

  fx <- ua_fixture(withr::local_tempdir())
  writeLines(c("##gff-version 3", "other_contig\tsrc\tCDS\t1\t90\t.\t+\t0\tName=COI"),
             file.path(fx$ann, "circ.gff"))
  expect_error(ua_init(fx), "seqid.*UA_CIRCULAR")

  fx <- ua_fixture(withr::local_tempdir(),
                   lin_gff = "CDS\t17900\t18100\t.\t+\t0\tID=c1;Name=COI")
  expect_error(ua_init(fx), "UA_LINEAR.*linear")

  fx <- ua_fixture(withr::local_tempdir(),
                   lin_gff = c("region\t1\t18024\t.\t+\t.\tIs_circular=true",
                               "CDS\t100\t300\t.\t+\t0\tName=COI"))
  expect_error(ua_init(fx), "Is_circular")

  fx <- ua_fixture(withr::local_tempdir(),
                   lin_gff = "CDS\t300\t100\t.\t+\t0\tID=c1;Name=COI")
  expect_error(ua_init(fx), "start after end|past the end|exceeds")

  fx <- ua_fixture(withr::local_tempdir())
  writeLines("mito_contig\t17900\t18100\tCOI\t0\t+", file.path(fx$ann, "lin.bed"))
  m <- utils::read.csv(fx$map)
  m$Annotation[1] <- "lin.bed"
  utils::write.csv(m, fx$map, row.names = FALSE)
  expect_error(ua_init(fx), "past the end|exceeds")

  fx <- ua_fixture(withr::local_tempdir(),
                   lin_gff = "repeat_unit\t100\t300\t.\t+\t.\tID=x")
  expect_error(ua_init(fx), "no recognized features")

  fx <- ua_fixture(withr::local_tempdir())
  m <- utils::read.csv(fx$map)
  m$Annotation[1] <- "lin.txt"
  m$Annotation[2] <- "missing.gff"
  utils::write.csv(m, fx$map, row.names = FALSE)
  err <- tryCatch(ua_init(fx), error = conditionMessage)
  expect_match(err, "lin.txt")
  expect_match(err, "missing.gff")

  fx <- ua_fixture(withr::local_tempdir())
  m <- utils::read.csv(fx$map)
  m$Annotation <- NULL
  utils::write.csv(m, fx$map, row.names = FALSE)
  expect_error(ua_init(fx), "'Annotation' column")
})

test_that("new_db_userAnnot rejects a multi-contig FASTA", {
  td <- withr::local_tempdir()
  fx <- ua_fixture(td, multi = TRUE)
  db <- file.path(td, ".sqlite")
  expect_error(suppressMessages(new_db_userAnnot(
    db_path = db, mapping_fn = fx$map, assembly_path = fx$asm,
    annotation_path = fx$ann, no_raw_data = TRUE)), "UA_MULTI_ONE.fasta.*one sequence")
  expect_false(file.exists(db))
})

test_that("force re-init removes stale per-sample annotation CSVs", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  fx <- ua_fixture(withr::local_tempdir())
  ua_init(fx)
  ua_dir <- file.path(fx$proj, "user_annotations")
  stale <- file.path(ua_dir, "UA_LINEAR_annotations.csv")
  writeLines("stale", stale)
  ua_init(fx, force = TRUE)
  expect_false(file.exists(stale))
  expect_true(dir.exists(file.path(ua_dir, "raw")))
  expect_length(list.files(file.path(ua_dir, "raw")), 2L)
})

test_that("reads are optional: no_raw_data = FALSE requires R1/R2", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  fx <- ua_fixture(withr::local_tempdir())
  expect_error(suppressMessages(new_project_userAnnot(
    path = fx$proj, mapping_fn = fx$map, assembly_path = fx$asm,
    annotation_path = fx$ann, data_path = fx$asm, Rproj = FALSE)), "R1")
})

test_that("unnamed features fall back to their note as raw name", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  fx <- ua_fixture(withr::local_tempdir(), lin_gff = c(
    "CDS\t100\t300\t.\t+\t0\tID=c1;Name=COI",
    "misc_feature\t4000\t4100\t.\t+\t.\tNote=control region",
    "misc_feature\t4200\t4300\t.\t+\t.\tNote=tandem repeat"
  ))
  ua_init(fx)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(fx$proj, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  gm <- DBI::dbGetQuery(con, "SELECT raw_name, gene FROM gene_name_map
                              WHERE feature_type = 'misc_feature' ORDER BY raw_name")
  expect_equal(gm$raw_name, c("control region", "tandem repeat"))
  expect_equal(gm$gene, c("ctrl", NA))
})

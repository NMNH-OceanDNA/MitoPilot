ua_resolve_project <- function(td) {
  src <- app_sys("test_data/assemblies")
  asm <- file.path(td, "asm")
  ann <- file.path(td, "ann")
  dir.create(asm)
  dir.create(ann)
  file.copy(file.path(src, c("UA_LINEAR.fasta", "UA_CIRCULAR.fasta")), asm)
  gff <- function(rows, fn) {
    writeLines(c("##gff-version 3", paste("mito_contig", "src", rows, sep = "\t")),
               file.path(ann, fn))
  }
  gff(c(
    "CDS\t100\t300\t.\t+\t0\tID=c1;Name=COX1",
    "tRNA\t700\t770\t.\t-\t.\tID=t1;Name=tRNA-Phe(gaa)",
    "CDS\t1000\t1200\t.\t+\t0\tID=o1;Name=ORF1",
    "CDS\t2000\t2100\t.\t+\t0\tID=n3;Name=ND3;ribosomal_slippage;Note=-1 frameshift",
    "CDS\t2100\t2400\t.\t+\t0\tID=n3;Name=ND3;ribosomal_slippage",
    "CDS\t3000\t3900\t.\t+\t0\tID=n2;Name=ND2;pseudo=true",
    "misc_feature\t4000\t4100\t.\t+\t.\tID=d1;Name=D-loop",
    "rRNA\t5000\t6000\t.\t+\t.\tID=r1;Name=16S"
  ), "lin.gff3")
  gff(c(
    "region\t1\t17265\t.\t+\t.\tID=reg;Is_circular=true",
    "CDS\t17000\t17400\t.\t+\t0\tID=c1;Name=cox1",
    "tRNA\t500\t570\t.\t+\t.\tID=t1;Name=trnL2(taa)"
  ), "circ.gff")
  m <- data.frame(ID = c("UA_LINEAR", "UA_CIRCULAR"), Taxon = "Danio rerio",
                  Assembly = c("UA_LINEAR.fasta", "UA_CIRCULAR.fasta"),
                  Topology = c("linear", "circular"), Annotation = c("lin.gff3", "circ.gff"))
  utils::write.csv(m, file.path(td, "mapping.csv"), row.names = FALSE)
  proj <- file.path(td, "proj")
  suppressMessages(new_project_userAnnot(
    path = proj, mapping_fn = file.path(td, "mapping.csv"), assembly_path = asm,
    annotation_path = ann, no_raw_data = TRUE, Rproj = FALSE, executor = "local"))
  proj
}

test_that("resolve, write, and confirm imported annotations", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_resolve_project(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  gm <- DBI::dbReadTable(con, "gene_name_map")
  expect_equal(sum(gm$raw_name %in% c("COX1", "cox1")), 2L)

  r <- resolve_user_annotations(con, "UA_LINEAR")
  expect_named(r, c("contig", "type", "gene", "product", "pos1", "pos2", "length",
                    "direction", "anticodon", "partial_start", "partial_stop",
                    "codon_start", "notes", "tool"))
  expect_true("cox1" %in% r$gene)
  expect_equal(r$product[r$gene == "cox1"], "cytochrome c oxidase subunit 1")
  expect_equal(r$product[r$gene == "trnF"], "tRNA-Phe")
  expect_equal(r$product[r$gene == "rrnL"], "16S ribosomal RNA")
  expect_false("ORF1" %in% r$gene)
  expect_equal(nrow(r), 6L)
  expect_true(all(r$partial_start %in% 0:1))
  expect_true(all(r$tool == "user"))
  expect_true(all(r$contig == "UA_LINEAR.1.1"))
  expect_false(is.unsorted(r$pos1))
  n3 <- r$notes[r$gene == "nad3"]
  expect_length(n3, 2L)
  expect_equal(n3, rep("JOIN: mode=frameshift group=1 note=-1 frameshift", 2))
  expect_match(r$notes[r$gene == "nad2"], "import: pseudogene")
  expect_equal(r$length[r$gene == "cox1"], 201L)
  expect_length(attr(r, "warnings"), 0L)

  rc <- resolve_user_annotations(con, "UA_CIRCULAR")
  w <- rc[rc$gene == "cox1", ]
  expect_true(w$pos1 > w$pos2)
  expect_equal(w$length, 401L)
  expect_equal(rc$product[rc$gene == "trnL"], "tRNA-Leu")

  changed <- write_user_annotation_csvs(con, td)
  expect_setequal(changed, c("UA_LINEAR", "UA_CIRCULAR"))
  fn <- file.path(td, "user_annotations", "UA_LINEAR_annotations.csv")
  csv <- readr::read_csv(fn, show_col_types = FALSE)
  expect_equal(nrow(csv), 6L)
  expect_false(any(grepl("(^|,)NA(,|$)", readLines(fn))))
  expect_length(write_user_annotation_csvs(con, td), 0L)

  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 0")
  set_gene_name_map(con, "ORF1", "CDS", "nad6")
  expect_false(user_annot_confirmed(con))
  out <- confirm_user_annotations(con, td)
  expect_true(user_annot_confirmed(con))
  expect_equal(out$requeued, "UA_LINEAR")
  expect_equal(out$kept, 9L)
  expect_equal(out$dropped, 1L)
  sw <- DBI::dbGetQuery(con, "SELECT ID, annotate_switch FROM annotate")
  expect_equal(sw$annotate_switch[sw$ID == "UA_LINEAR"], 1L)
  expect_equal(sw$annotate_switch[sw$ID == "UA_CIRCULAR"], 0L)
  expect_false(any(is.na(DBI::dbReadTable(con, "gene_name_map")$time_stamp)))

  set_gene_name_map(con, "COX1", "CDS", "(drop)")
  gm <- DBI::dbReadTable(con, "gene_name_map")
  expect_true(is.na(gm$gene[gm$raw_name == "COX1"]))
  expect_equal(gm$manual[gm$raw_name == "COX1"], 1L)
  expect_false(user_annot_confirmed(con))
  reset_gene_name_map(con)
  gm <- DBI::dbReadTable(con, "gene_name_map")
  expect_equal(gm$gene[gm$raw_name == "COX1"], "cox1")
  expect_true(is.na(gm$gene[gm$raw_name == "ORF1"]))
  expect_true(all(gm$manual == 0L & gm$confirmed == 0L))
})

test_that("split JOIN mapping drops the marker with a warning; empty sample gets header", {
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible(NULL))
  td <- ua_resolve_project(withr::local_tempdir())
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE user_annotations SET raw_name = 'ND3x'
                       WHERE ID = 'UA_LINEAR' AND raw_name = 'ND3' AND pos1 = 2100")
  DBI::dbExecute(con, "INSERT INTO gene_name_map (raw_name, feature_type, gene, guess)
                       VALUES ('ND3x', 'CDS', 'nad6', 'nad6')")
  r <- resolve_user_annotations(con, "UA_LINEAR")
  expect_false(any(grepl("JOIN:", r$notes)))
  expect_equal(r$notes[r$gene == "nad3"], "import: -1 frameshift")
  expect_match(attr(r, "warnings"), "ND3")

  DBI::dbExecute(con, "UPDATE gene_name_map SET gene = NULL")
  expect_equal(nrow(resolve_user_annotations(con, "UA_CIRCULAR")), 0L)
  write_user_annotation_csvs(con, td)
  txt <- readLines(file.path(td, "user_annotations", "UA_CIRCULAR_annotations.csv"))
  expect_equal(txt, "contig,type,gene,product,pos1,pos2,length,direction,anticodon,partial_start,partial_stop,codon_start,notes,tool")
})

test_that("NCBI GFF3 collapses hierarchy and prefers gene= over Name=", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/ncbi.gff3"))
  nd6 <- x[x$raw_name == "ND6", ]
  expect_equal(c(nd6$partial_start, nd6$partial_stop), c(0L, 1L))   # minus strand: low end open = 3'
  expect_equal(x$partial_start[x$raw_name == "COX1"], 0L)
  expect_setequal(x$feature_type, c("CDS", "tRNA", "rRNA", "D_loop"))
  expect_true("COX1" %in% x$raw_name)
  expect_false(any(grepl("^YP_", x$raw_name)))
  expect_true(isTRUE(attr(x, "is_circular")))
  expect_equal(unname(attr(x, "dropped_types")[["exon"]]), 1L)
})
test_that("messy GFF3 merges two-part origin-wrapping CDS into pos1 > pos2", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/messy.gff3"), contig_length = 16500)
  cob <- x[x$feature_id == "cds-cob", ]
  expect_equal(nrow(cob), 1L); expect_gt(cob$pos1, cob$pos2)
  expect_equal(sum(x$raw_name == "ND1"), 1L)   # gene+CDS same span -> one row typed CDS
  expect_equal(x$partial_start[grepl("^COI", x$raw_name)], 1L)   # "5' Partial" in Name
  nad5 <- x[x$feature_id == "cds-nad5", ]
  expect_equal(nrow(nad5), 2L); expect_equal(unique(nad5$segment_mode), "exon"); expect_length(unique(nad5$segment_group), 1L)
  nad3 <- x[x$feature_id == "cds-nad3", ]
  expect_equal(unique(nad3$segment_mode), "frameshift"); expect_equal(unique(nad3$segment_note), "-1 frameshift")
  expect_true(is.na(cob$segment_group))                     # wrap pair is one segment, not a group
  expect_equal(x$pseudo[x$raw_name == "ND2"], 1L); expect_equal(x$note[x$raw_name == "ND2"], "frameshifted")
  expect_equal(x$codon_start[x$raw_name == "ND4"], 2L); expect_equal(x$partial_start[x$raw_name == "ND4"], 1L)
  expect_true(jsonlite::fromJSON(x$attrs[x$raw_name == "ATP8"])$transl_except != "")
  expect_false(any(grepl("transl_except", attr(x, "warnings"))))   # length not multiple of 3: consistent
})
test_that("BED12 blocks become exon segments", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/blocks.bed"))
  expect_equal(nrow(x[x$raw_name == "ND5", ]), 2L); expect_equal(unique(x$segment_mode[x$raw_name == "ND5"]), "exon")
  expect_equal(x$pos1[x$raw_name == "ND5"][1], 101L)   # chromStart 100 + blockStart 0 + 1
})
test_that("BED is converted to 1-based and defaults strand to +", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/simple.bed"))
  expect_equal(x$pos1[x$raw_name == "COX1"], 1L)
  expect_equal(x$pos2[x$raw_name == "COX1"], 1545L)
  expect_true(all(is.na(x$feature_type)))
  expect_equal(x$direction[3], "+")
})
test_that("CRLF and trailing tabs parse like LF", {
  a <- parse_user_annotation(test_path("fixtures/user_annot/mitos.gff"))
  b <- parse_user_annotation(test_path("fixtures/user_annot/crlf.gff3"))
  expect_equal(a$raw_name, b$raw_name); expect_equal(a$pos1, b$pos1)
})
test_that("MitoPilot export .tbl round-trips: strand, pairs, wrap, anticodon note", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/mitopilot_export.tbl"), contig_length = 16500)
  expect_false("gene" %in% x$feature_type)                     # gene rows collapsed into children
  expect_equal(x$direction[x$raw_name == "nad6"], "-")
  expect_lt(x$pos1[x$raw_name == "nad6"], x$pos2[x$raw_name == "nad6"])   # stored lo..hi
  cob <- x[x$raw_name == "cob", ]; expect_equal(c(cob$pos1, cob$pos2), c(16400L, 700L))
  expect_equal(x$anticodon_raw[x$raw_name == "trnF"], "anticodon:gaa")
  expect_equal(unname(attr(x, "dropped_types")[["gap"]]), 1L)
  expect_equal(c(x$partial_start[x$raw_name == "nad1"], x$partial_stop[x$raw_name == "nad1"]), c(1L, 0L))
  expect_equal(c(x$partial_start[x$raw_name == "rrn16"], x$partial_stop[x$raw_name == "rrn16"]), c(0L, 1L))
  expect_equal(x$pos1[x$raw_name == "nad1"], 2000L)   # marker stripped from the coordinate
  expect_equal(x$codon_start[x$raw_name == "nad1"], 2L)
  nad3 <- x[x$raw_name == "nad3", ]; expect_equal(nrow(nad3), 2L); expect_equal(unique(nad3$segment_mode), "frameshift")
  expect_equal(unique(nad3$segment_note), "-1 frameshift")
  expect_equal(x$pseudo[x$note %in% "frameshifted"], 1L)
})
test_that("NCBI .tbl: gb|ACC| seqid reduced, minus-strand 3' partial, gene > product", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/ncbi.tbl"))
  expect_equal(unique(x$seqid), "NC_000001.1")
  expect_true("COX1" %in% x$raw_name[x$feature_type == "CDS"])
  nd6 <- x[x$raw_name == "ND6", ]
  expect_equal(nd6$direction, "-"); expect_equal(c(nd6$partial_start, nd6$partial_stop), c(0L, 1L))
  expect_true(is.integer(x$pos1) && is.integer(x$pos2))
})
test_that("multiple seqids is an error", {
  tmp <- withr::local_tempfile(fileext = ".gff3")
  writeLines(c("chrA\t.\tCDS\t1\t9\t.\t+\t.\tName=x", "chrB\t.\tCDS\t1\t9\t.\t+\t.\tName=y"), tmp)
  expect_error(parse_user_annotation(tmp), "seqid")
})
test_that("GFF3: literal % decodes, end past contig length wraps, extra kept types", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/messy.gff3"), contig_length = 16500)
  expect_equal(x$note[x$feature_id %in% "rep1"], "50% identity")
  expect_equal(c(x$pos1[x$raw_name == "trnP"], x$pos2[x$raw_name == "trnP"]), c(16450L, 20L))
})
test_that(".tbl: exported JOIN gene span collapses with partials; minus wrap keeps strand", {
  x <- parse_user_annotation(test_path("fixtures/user_annot/mitopilot_export.tbl"), contig_length = 16500)
  nad3 <- x[x$raw_name %in% "nad3", ]
  expect_equal(nrow(nad3), 2L)
  expect_equal(nad3$feature_type, c("CDS", "CDS"))
  expect_equal(nad3$partial_stop, c(0L, 1L))
  expect_equal(nad3$partial_start, c(0L, 0L))
  nad5 <- x[x$raw_name %in% "nad5", ]
  expect_equal(nrow(nad5), 1L)
  expect_equal(nad5$direction, "-")
  expect_equal(c(nad5$pos1, nad5$pos2), c(16480L, 1L))
  tmp <- withr::local_tempfile(fileext = ".tbl")
  writeLines(c(">Feature x", "1\t50\trepeat_region", "60\t70\torigin_of_replication"), tmp)
  expect_equal(nrow(parse_user_annotation(tmp)), 2L)
})

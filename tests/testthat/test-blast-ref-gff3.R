# tests/testthat/test-blast-ref-gff3.R
gff3_lines <- function(...) paste(c(
  "##gff-version 3", "##sequence-region NC_1 1 16500",
  "NC_1\tGenbank\tregion\t1\t16500\t.\t+\t.\tID=NC_1:1..16500;Is_circular=true;gbkey=Src",
  ...), collapse = "\n")

test_that("reference GFF3 partial ranges become 5'/3' partial flags by strand", {
  g <- .parse_ref_gff3("NC_1", gff3_lines(
    "NC_1\tGenbank\tCDS\t1\t900\t.\t+\t.\tID=cds1;gene=ND1;gbkey=CDS;partial=true;start_range=.,1;transl_table=2",
    "NC_1\tGenbank\tCDS\t1000\t2100\t.\t-\t.\tID=cds2;gene=COX1;gbkey=CDS;partial=true;start_range=.,1000",
    "NC_1\tGenbank\tCDS\t3000\t4100\t.\t+\t.\tID=cds3;gene=CYTB;gbkey=CDS;partial=true;end_range=4100,.",
    "NC_1\tGenbank\ttRNA\t5000\t5070\t.\t+\t.\tID=rna1;product=tRNA-Phe;gbkey=tRNA"))
  r <- g$result
  expect_equal(names(r), c("gene", "type", "pos1", "pos2", "direction", "ref_length",
                           "partial_start", "partial_stop"))
  expect_equal(r$partial_start, c(1L, 0L, 0L, 0L))
  expect_equal(r$partial_stop, c(0L, 1L, 1L, 0L))
})

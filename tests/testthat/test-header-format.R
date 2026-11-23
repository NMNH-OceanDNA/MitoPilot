d <- data.frame(ID = c("s1", "s2"), seqid = c("s1", "s2"), Taxon = "Gadus morhua",
                loc = c("Norway", "Faroe [north]"), stringsAsFactors = FALSE)
lvl <- function(t, data = d[1, ]) validate_fasta_header(t, data)$level
msg <- function(t, data = d[1, ]) validate_fasta_header(t, data)$message

test_that("header_clean gives one regularly spaced line", {
  expect_equal(
    header_clean("s1\n[ organism = Gadus morhua ][note=a\tb]\r\n  title  "),
    "s1 [organism=Gadus morhua] [note=a b] title"
  )
  expect_equal(header_fill(d[1, ], "{seqid} [organism={Taxon}]\n{Taxon}"),
               "s1 [organism=Gadus morhua] Gadus morhua")
})

test_that("a line break in the template is not an error", {
  expect_equal(lvl("{seqid} [organism={Taxon}]\n{Taxon}"), "ok")
})

test_that("broken brackets block", {
  expect_equal(lvl("{seqid} [organism={Taxon} {Taxon}"), "error")
  expect_equal(lvl("{seqid} organism={Taxon}] {Taxon}"), "error")
  expect_equal(lvl("{seqid} [organism] {Taxon}"), "error")
  expect_match(msg("{seqid} [organism={Taxon}] [geo_loc_name={loc}]", d), "s2")
  expect_equal(lvl("{seqid} [organism={Taxon}] [geo_loc_name={loc}]", d), "error")
})

test_that("format slips warn", {
  expect_match(msg("x_{seqid} [organism={Taxon}]"), "Start the header")
  expect_match(msg("{seqid} [organism={Taxon}] [country=Norway]"), "geo-loc-name")
  expect_match(msg("{seqid} [organism={Taxon}] [Voucher URI=x]"), "not on NCBI")
  expect_match(msg("{seqid} [organism={Taxon}] [authority=L.]"), "no longer accepts")
  expect_match(msg("{seqid} {Taxon} [organism={Taxon}]"), "title text")
  expect_match(msg("{seqid} [organism={Taxon}] mito [note=x]"), "title text")
  long <- transform(d[1, ], seqid = strrep("a", 26))
  expect_match(msg("{seqid} [organism={Taxon}]", long), "25 characters")
  expect_match(msg("{seqid} [organism={Taxon}]", transform(d[1, ], seqid = "a/b")), "25 characters")
  expect_match(msg("{seqid} [organism={Taxon}]", transform(d[1, ], Taxon = "Gadus mörhua")),
               "non-ASCII")
  expect_equal(lvl("x_{seqid} [organism={Taxon}]"), "warn")
  expect_match(msg("{seqid}[organism={Taxon}]"), "Start the header")
  expect_equal(lvl("{seqid}_x [organism={Taxon}]"), "ok")
})

test_that("accepted spellings pass", {
  expect_equal(lvl("{seqid} [Organism={Taxon}] [specimen_voucher=a] [lat-lon=1 N 2 W] [BioSample=x] t"), "ok")
})

test_that("a token that names an R function reads as an unknown field", {
  expect_match(msg("{seq} [organism={Taxon}]"), "\\{seq\\} is not one of the available fields")
})

test_that("header_lock puts the record name first", {
  expect_equal(header_lock("{seqid} [organism=x] t"), "{seqid} [organism=x] t")
  expect_equal(header_lock("{ID} [organism=x] t"), "{seqid} [organism=x] t")
  expect_equal(header_lock("[organism=x] t"), "{seqid} [organism=x] t")
  expect_equal(header_lock("{seqid} [organism=x]", "{seqid}_{element}"), "{seqid}_{element} [organism=x]")
})

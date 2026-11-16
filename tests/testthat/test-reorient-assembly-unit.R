# "Rotate" in the Annotate details window: rotate a circular unit to start on
# a chosen feature, flipping the strand when that feature is on the minus strand.

reorient_db <- function(topology = "circular", env = parent.frame()) {
  set.seed(7)
  seq_chr <- paste(sample(c("A", "C", "G", "T"), 60, TRUE), collapse = "")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  withr::defer(DBI::dbDisconnect(con), envir = env)
  DBI::dbExecute(con, "CREATE TABLE assemblies (
    ID TEXT, path INTEGER, scaffold INTEGER, sequence TEXT, length INTEGER, topology TEXT,
    depth TEXT, gc TEXT, errors TEXT, edited INTEGER, time_stamp INTEGER)")
  DBI::dbExecute(con, "CREATE TABLE annotations (
    ID TEXT, path INTEGER, scaffold INTEGER, gene TEXT, pos1 INTEGER, pos2 INTEGER, direction TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate (
    ID TEXT, path INTEGER, scaffold INTEGER, length INTEGER, topology TEXT, structure TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES ('S1',1,1,?,60,?,?,NULL,NULL,0,0)",
                 params = list(seq_chr, topology, paste(1:60, collapse = " ")))
  DBI::dbExecute(con, "INSERT INTO annotate VALUES ('S1',1,1,60,?,'geneA|geneB|geneC')",
                 params = list(topology))
  DBI::dbAppendTable(con, "annotations", data.frame(
    ID = "S1", path = 1L, scaffold = 1L,
    gene = c("geneA", "geneB", "geneC", "gone"),
    pos1 = c(5L, 21L, 55L, 0L), pos2 = c(14L, 40L, 3L, 0L),
    direction = c("+", "-", "+", "+")
  ))
  list(con = con, seq = seq_chr, dir = withr::local_tempdir(.local_envir = env))
}

feat_seq <- function(seq_chr, p1, p2, direction) {
  s <- if (p1 <= p2) substr(seq_chr, p1, p2) else paste0(substr(seq_chr, p1, nchar(seq_chr)), substr(seq_chr, 1, p2))
  if (direction == "-") s <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))
  s
}

# every live feature must still read the same bases, 5' to 3'
expect_features_kept <- function(db) {
  before <- data.frame(gene = c("geneA", "geneB", "geneC"),
                       pos1 = c(5L, 21L, 55L), pos2 = c(14L, 40L, 3L), direction = c("+", "-", "+"))
  now <- DBI::dbGetQuery(db$con, "SELECT * FROM annotations WHERE pos1 > 0")
  new_seq <- DBI::dbGetQuery(db$con, "SELECT sequence FROM assemblies")$sequence
  for (i in seq_len(nrow(before))) {
    j <- match(before$gene[i], now$gene)
    expect_equal(
      feat_seq(new_seq, now$pos1[j], now$pos2[j], now$direction[j]),
      feat_seq(db$seq, before$pos1[i], before$pos2[i], before$direction[i])
    )
  }
  now
}

test_that("a plus-strand feature becomes position 1 by rotation alone", {
  db <- reorient_db()
  res <- reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "geneA", 5)
  expect_equal(res$rotated, 4L)
  expect_false(res$flipped)
  now <- expect_features_kept(db)
  expect_equal(now$pos1[now$gene == "geneA"], 1L)
  expect_equal(now$pos1[now$gene == "geneC"], 51L)  # still spans the origin
  expect_equal(now$pos2[now$gene == "geneC"], 59L)
  asm <- DBI::dbGetQuery(db$con, "SELECT * FROM assemblies")
  expect_equal(asm$depth, paste(c(5:60, 1:4), collapse = " "))
  expect_equal(asm$edited, 1L)
  expect_equal(DBI::dbGetQuery(db$con, "SELECT structure FROM annotate")$structure,
               "geneA|geneB|geneC")
  # tombstone untouched
  expect_equal(DBI::dbGetQuery(db$con, "SELECT pos1, pos2 FROM annotations WHERE gene = 'gone'"),
               data.frame(pos1 = 0L, pos2 = 0L))
})

test_that("a minus-strand feature flips the unit and ends up first on the plus strand", {
  db <- reorient_db()
  res <- reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "geneB", 21)
  expect_true(res$flipped)
  now <- expect_features_kept(db)
  b <- now[now$gene == "geneB", ]
  expect_equal(c(b$pos1, b$pos2, b$direction), c("1", "20", "+"))
  expect_equal(now$direction[now$gene == "geneA"], "-")
  # old positions 40..1, then 60..41
  expect_equal(DBI::dbGetQuery(db$con, "SELECT depth FROM assemblies")$depth,
               paste(c(40:1, 60:41), collapse = " "))
  expect_equal(DBI::dbGetQuery(db$con, "SELECT structure FROM annotate")$structure,
               "geneB|geneA|geneC")
})

test_that("restore undoes a reorient", {
  db <- reorient_db()
  reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "geneB", 21)
  expect_equal(assembly_backup_ops(db$con, "S1", 1, 1), "reorient")
  restore_assembly_unit(db$con, "S1", 1, 1, db$dir)
  expect_equal(DBI::dbGetQuery(db$con, "SELECT sequence FROM assemblies")$sequence, db$seq)
  expect_features_kept(db)
})

test_that("reorient refuses linear units, no-ops and cuts through a feature", {
  lin <- reorient_db("linear")
  expect_error(reorient_assembly_unit(lin$con, "S1", 1, 1, lin$dir, "geneA", 5), "circular")
  db <- reorient_db()
  expect_error(reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "nope", 5), "not found")
  DBI::dbExecute(db$con, "INSERT INTO annotations VALUES ('S1',1,1,'over',4,8,'+')")
  expect_error(reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "geneA", 5), "split the over")
  expect_equal(DBI::dbGetQuery(db$con, "SELECT sequence FROM assemblies")$sequence, db$seq)
  reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "over", 4)
  expect_error(reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "over", 1), "already starts")
})

test_that("the coverage file and per-base tracks follow a flip", {
  db <- reorient_db()
  cov_fn <- unit_annotate_file(db$dir, "S1", 1, 1, "coverageStats")
  dir.create(dirname(cov_fn), recursive = TRUE)
  utils::write.csv(data.frame(
    SeqId = "S1.1.1", Position = 1:60, Call = strsplit(db$seq, "")[[1]],
    MeanDepth = 101:160, GC = 0.5, ErrorRate = 0
  ), cov_fn, row.names = FALSE, quote = FALSE)
  reorient_assembly_unit(db$con, "S1", 1, 1, db$dir, "geneB", 21)
  cov <- utils::read.csv(cov_fn)
  asm <- DBI::dbGetQuery(db$con, "SELECT sequence, depth FROM assemblies")
  expect_equal(cov$Position, 1:60)
  expect_equal(paste(cov$Call, collapse = ""), asm$sequence)
  expect_equal(cov$MeanDepth, c(140:101, 160:141))
  expect_equal(asm$depth, paste(c(140:101, 160:141), collapse = " "))
  expect_match(readLines(unit_annotate_file(db$dir, "S1", 1, 1, "assembly"))[1], "circular")
})

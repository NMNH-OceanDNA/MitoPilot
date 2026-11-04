# "Found by MitoPilot": features added by the annotation tools in a
# user-annotation project, counted per unit for the Annotate table.

tc_con <- function(env = parent.frame()) {
  td <- withr::local_tempdir(.local_envir = env)
  mapping <- file.path(td, "mapping.csv")
  utils::write.csv(
    data.frame(ID = c("U", "V"), Taxon = "Danio rerio", R1 = "a", R2 = "b"),
    mapping, row.names = FALSE
  )
  suppressMessages(new_db(db_path = file.path(td, ".sqlite"), mapping_fn = mapping))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(td, ".sqlite"))
  withr::defer(DBI::dbDisconnect(con), envir = env)
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  DBI::dbExecute(con, "UPDATE annotate SET path = 1, scaffold = 1")
  DBI::dbAppendTable(con, "assemblies", data.frame(
    ID = c("U", "V"), path = 1L, scaffold = 1L, length_raw = 16000L, ignore = 0L
  ))
  con
}

tc_rows <- function(con, tool, gene, type, pos1) {
  DBI::dbAppendTable(con, "annotations", data.frame(
    ID = "U", path = 1L, scaffold = 1L, tool = tool, gene = gene,
    type = type, pos1 = as.integer(pos1)
  ))
}

test_that("toolCount counts non-user, non-ORF rows per unit", {
  con <- tc_con()
  tc_rows(con,
    tool = c("user", "user", "user", "tRNAscan-SE", "MITOS2", "ORFfinder"),
    gene = c("cox1", "rrnL", "nad1", "trnC", "trnA", "ORF1"),
    type = c("PCG", "rRNA", "PCG", "tRNA", "tRNA", "ORF"),
    pos1 = c(10, 200, 900, 700, 500, 950)
  )
  d <- fetch_annotate_data(session = list(userData = list(con = con, user_annot = TRUE)))
  expect_identical(d$toolCount[d$ID == "U"], 2L)
  expect_equal(d$toolGenes[d$ID == "U"], "trnA; trnC")
  expect_equal(which(names(d) == "toolCount"), which(names(d) == "ORFCount") + 1L)
})

test_that("toolCount is 0 (not NA) with no annotations and in a tools-only project", {
  con <- tc_con()
  tc_rows(con, tool = "MITOS2", gene = "trnA", type = "tRNA", pos1 = 500)
  d <- fetch_annotate_data(session = list(userData = list(con = con, user_annot = TRUE)))
  expect_identical(d$toolCount[d$ID == "V"], 0L)
  expect_equal(d$toolGenes[d$ID == "V"], "")
  d <- fetch_annotate_data(session = list(userData = list(con = con, user_annot = FALSE)))
  expect_identical(d$toolCount, c(0L, 0L))
  expect_equal(d$toolGenes, c("", ""))
})

test_that("toolCount is in the Counts group after ORFCount", {
  g <- ANNOTATE_COL_GROUPS$Counts
  expect_equal(which(g == "toolCount"), which(g == "ORFCount") + 1L)
  expect_equal(MP_COL_NAMES[["toolCount"]], "Found by MitoPilot")
})

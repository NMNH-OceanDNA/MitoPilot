# tests/testthat/test-user-annot-track.R
rot <- function(s, n) paste0(substr(s, n + 1, nchar(s)), substr(s, 1, n))
set.seed(1)
S <- paste(sample(c("A", "C", "G", "T"), 2000, TRUE), collapse = "")

test_that("find_frame: identity, rotation, rc, probe over the new origin", {
  expect_equal(find_frame(S, S), list(offset = 0L, rc = FALSE))
  expect_equal(find_frame(S, rot(S, 500)), list(offset = -500L, rc = FALSE))
  expect_equal(find_frame(S, revcomp_chr(S)), list(offset = 0L, rc = TRUE))
  expect_equal(find_frame(S, rot(S, 1990)), list(offset = -1990L, rc = FALSE))
  expect_equal(find_frame(S, rot(S, 30)), list(offset = -30L, rc = FALSE))
  expect_equal(find_frame(S, revcomp_chr(rot(S, 500))), list(offset = -500L, rc = TRUE))
})

test_that("find_frame: trimmed start uses a later probe; repeats and unrelated give NULL", {
  expect_equal(find_frame(S, substr(S, 101, 1900)), list(offset = -100L, rc = FALSE))
  expect_null(find_frame(S, paste0(S, S)))
  expect_null(find_frame(S, paste(rev(strsplit(S, "")[[1]]), collapse = "")))
})

test_that("map_feature: rotation keeps strand, rc flips, wrap allowed", {
  f <- list(offset = -500L, rc = FALSE)
  expect_equal(map_feature(600L, 700L, "+", f, 2000L, 2000L)[c("pos1", "pos2", "dir")],
               list(pos1 = 100L, pos2 = 200L, dir = "+"))
  expect_equal(map_feature(400L, 600L, "+", f, 2000L, 2000L)[c("pos1", "pos2")],
               list(pos1 = 1900L, pos2 = 100L))
  r <- list(offset = 0L, rc = TRUE)
  expect_equal(map_feature(10L, 30L, "+", r, 2000L, 2000L)[c("pos1", "pos2", "dir")],
               list(pos1 = 1971L, pos2 = 1991L, dir = "-"))
  expect_false(map_feature(10L, 30L, "+", r, 2000L, 2000L)$clipped)
})

test_that("map_feature: linear shift clips, flags, and drops", {
  f <- list(offset = -100L, rc = FALSE)
  expect_equal(map_feature(50L, 150L, "+", f, 2000L, 1800L),
               list(pos1 = 1L, pos2 = 50L, dir = "+", clipped = TRUE))
  expect_false(map_feature(500L, 600L, "+", f, 2000L, 1800L)$clipped)
  expect_null(map_feature(10L, 90L, "+", f, 2000L, 1800L))
  expect_equal(map_feature(1950L, 300L, "+", f, 2000L, 1800L)[c("pos1", "pos2", "clipped")],
               list(pos1 = 1L, pos2 = 200L, clipped = TRUE))
})

track_db <- function(orig, cur, backup = NULL, ua, gm, with_contigs = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies (ID TEXT, path INTEGER, scaffold INTEGER,
                       sequence TEXT, length INTEGER, topology TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES ('S1', 1, 1, ?, ?, 'circular')",
                 params = list(cur, nchar(cur)))
  if (with_contigs) {
    DBI::dbExecute(con, "CREATE TABLE user_contigs (ID TEXT PRIMARY KEY, sequence TEXT, length INTEGER)")
    DBI::dbExecute(con, "INSERT INTO user_contigs VALUES ('S1', ?, ?)", params = list(orig, nchar(orig)))
  }
  if (!is.null(backup)) {
    ensure_assembly_backup(con)
    DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length, ops)
                         VALUES ('S1', 1, 1, ?, ?, 'linearize+trim')",
                   params = list(backup, nchar(backup)))
  }
  DBI::dbWriteTable(con, "user_annotations", ua)
  DBI::dbWriteTable(con, "gene_name_map", gm)
  con
}
ua_rows <- function(raw_name, pos1, pos2, direction = "+", feature_type = "CDS") {
  data.frame(ID = "S1", feature_idx = seq_along(raw_name), raw_name = raw_name,
             feature_type = feature_type, pos1 = as.integer(pos1), pos2 = as.integer(pos2),
             direction = direction, segment_group = NA_integer_, stringsAsFactors = FALSE)
}
gm_rows <- function(raw_name, gene, feature_type = "CDS") {
  data.frame(raw_name = raw_name, feature_type = feature_type, gene = gene, stringsAsFactors = FALSE)
}

test_that("original_features: no edits, WF2 rotation only", {
  con <- track_db(S, rot(S, 500), ua = ua_rows("COI", 600, 700), gm = gm_rows("COI", "cox1"))
  withr::defer(DBI::dbDisconnect(con))
  o <- original_features(con, "S1", 1L, 1L)
  expect_true(o$ok)
  expect_equal(o$length, 2000L)
  r <- o$rows
  expect_equal(r$name, "COI"); expect_equal(r$gene, "cox1"); expect_equal(r$type, "PCG")
  expect_equal(c(r$pos1, r$pos2), c(100L, 200L)); expect_equal(r$dir, "+")
  expect_equal(c(r$orig1, r$orig2), c(600L, 700L)); expect_equal(r$odir, "+")
  expect_false(r$clipped)
})

test_that("original_features: backup after Linearize, then Trim", {
  mid <- rot(S, 500)
  con <- track_db(S, substr(mid, 51, 1800), backup = mid,
                  ua = ua_rows(c("COI", "ND1"), c(600, 520), c(700, 580)),
                  gm = gm_rows(c("COI", "ND1"), c("cox1", "nad1")))
  withr::defer(DBI::dbDisconnect(con))
  o <- original_features(con, "S1", 1L, 1L)
  expect_true(o$ok)
  expect_equal(o$length, 1750L)
  r <- o$rows
  expect_equal(r$pos1[r$name == "COI"], 50L); expect_equal(r$pos2[r$name == "COI"], 150L)
  expect_false(r$clipped[r$name == "COI"])
  expect_equal(r$pos1[r$name == "ND1"], 1L); expect_equal(r$pos2[r$name == "ND1"], 30L)
  expect_true(r$clipped[r$name == "ND1"])
})

test_that("original_features: unmapped names keep their raw type; missing table or frame give reasons", {
  con <- track_db(S, S, ua = ua_rows(c("x", "y", "z"), c(10, 50, 90), c(20, 60, 100),
                                     feature_type = c("tRNA", "misc_feature", "CDS")),
                  gm = gm_rows(c("x", "y", "z"), c(NA, "ctrl", NA), c("tRNA", "misc_feature", "CDS")))
  withr::defer(DBI::dbDisconnect(con))
  r <- original_features(con, "S1", 1L, 1L)$rows
  expect_equal(r$type, c("tRNA", "ctrl", "PCG"))
  expect_equal(is.na(r$gene), c(TRUE, FALSE, TRUE))

  con2 <- track_db(S, S, ua = ua_rows("COI", 1, 9), gm = gm_rows("COI", "cox1"), with_contigs = FALSE)
  withr::defer(DBI::dbDisconnect(con2))
  o <- original_features(con2, "S1", 1L, 1L)
  expect_false(o$ok); expect_equal(o$reason, "Original sequence not stored for this sample")

  con3 <- track_db(S, paste(rev(strsplit(S, "")[[1]]), collapse = ""),
                   ua = ua_rows("COI", 1, 9), gm = gm_rows("COI", "cox1"))
  withr::defer(DBI::dbDisconnect(con3))
  o <- original_features(con3, "S1", 1L, 1L)
  expect_false(o$ok)
  expect_equal(o$reason, "Original annotation cannot be placed on the edited sequence")
})

test_that("flag_original: one status each; name-only change is unchanged; duplicates pair one-to-one", {
  rows <- data.frame(
    name = c("COI", "ND1", "ND2", "REP", "ATP6", "trnL", "trnL", "CYTB"),
    gene = c("cox1", "nad1", "nad2", NA, "atp6", "trnL", "trnL", "cob"),
    pos1 = c(100L, 300L, 500L, 700L, 800L, 1000L, 1200L, 1500L),
    pos2 = c(200L, 403L, 600L, 750L, 900L, 1070L, 1270L, 1600L),
    dir = c("+", "+", "+", "+", "+", "+", "+", "-"),
    clipped = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  curated <- data.frame(
    gene = c("cox1", "nad1", "nad2", "atp6", "trnL", "trnL", "cob"),
    pos1 = c(100L, 300L, 0L, 800L, 1200L, 1000L, 1500L),
    pos2 = c(200L, 400L, 0L, 900L, 1275L, 1070L, 1600L),
    direction = c("+", "+", "+", "+", "+", "+", "-"),
    stringsAsFactors = FALSE
  )
  s <- flag_original(rows, curated, 2000L)$status
  expect_equal(s, c("unchanged", "moved", "dropped", "unmapped", "clipped",
                    "unchanged", "moved", "unchanged"))
})

test_that("flag_original: wrapped original pairs with its wrapped curated row", {
  rows <- data.frame(name = "COI", gene = "cox1", pos1 = 1950L, pos2 = 100L, dir = "+",
                     clipped = FALSE, stringsAsFactors = FALSE)
  curated <- data.frame(gene = "cox1", pos1 = 1950L, pos2 = 100L, direction = "+",
                        stringsAsFactors = FALSE)
  expect_equal(flag_original(rows, curated, 2000L)$status, "unchanged")
  curated$pos2 <- 103L
  expect_equal(flag_original(rows, curated, 2000L)$status, "moved")
})

test_that("viewer UIs render the Original checkbox only for user-annotation projects", {
  sv <- as.character(seqview_ui("x", user_annot = TRUE))
  expect_match(sv, "x-show_orig"); expect_match(sv, "x-orig_note")
  expect_no_match(as.character(seqview_ui("x")), "show_orig")
  sy <- as.character(synteny_viewer_ui(shiny::NS("y"), user_annot = TRUE))
  expect_match(sy, "y-show_orig"); expect_match(sy, "y-synteny_orig_note")
  expect_no_match(as.character(synteny_viewer_ui(shiny::NS("y"))), "show_orig")
})

test_that("map_feature: rotation then trim places features on both sides of the cut", {
  # Linearize rotated by 700 then Trim cut 200 bp off the end; probe at 1 gives +1300
  f <- list(offset = 1300L, rc = FALSE)
  expect_equal(map_feature(100L, 200L, "+", f, 2000L, 1800L)[c("pos1", "pos2", "clipped")],
               list(pos1 = 1400L, pos2 = 1500L, clipped = FALSE))
  expect_equal(map_feature(900L, 1000L, "+", f, 2000L, 1800L)[c("pos1", "pos2", "clipped")],
               list(pos1 = 200L, pos2 = 300L, clipped = FALSE))
  expect_null(map_feature(550L, 650L, "+", f, 2000L, 1800L))          # inside the trimmed 200 bp
  expect_equal(map_feature(450L, 560L, "+", f, 2000L, 1800L)[c("pos1", "pos2", "clipped")],
               list(pos1 = 1750L, pos2 = 1800L, clipped = TRUE))
})

test_that("original_features: backup is the pre-Linearize sequence; Linearize + Trim keeps every feature", {
  cur <- substr(rot(S, 700), 1, 1800)
  con <- track_db(S, cur, backup = S,
                  ua = ua_rows(c("A", "B", "C"), c(100, 900, 550), c(200, 1000, 650)),
                  gm = gm_rows(c("A", "B", "C"), c("cox1", "nad1", "nad2")))
  withr::defer(DBI::dbDisconnect(con))
  r <- original_features(con, "S1", 1L, 1L)$rows
  expect_equal(r$name, c("A", "B"))
  expect_equal(r$pos1, c(1400L, 200L)); expect_equal(r$pos2, c(1500L, 300L))
})

test_that("map_feature: an origin-spanning original becomes one contiguous feature after Linearize + Trim", {
  # GB_MITOS cob 15912..523 on a 16529 bp contig; Linearize before trnF (1530) then Trim to 15825
  f <- list(offset = 15000L, rc = FALSE)
  expect_equal(map_feature(15912L, 523L, "+", f, 16529L, 15825L),
               list(pos1 = 14383L, pos2 = 15523L, dir = "+", clipped = FALSE))
})

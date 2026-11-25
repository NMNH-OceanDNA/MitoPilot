codes <- function(f, ...) {
  r <- assembly_pass(f$con, ...)$reasons
  if (nzchar(r)) strsplit(r, ",")[[1]] else character()
}

test_that("clean circular single-path sample passes and is locked as auto", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80)
  expect_true(assembly_pass(f$con)$pass)
  r <- auto_accept(f$con, f$dir_out, "assemble")
  expect_equal(r$accepted$ID, "A")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT action FROM review_log")[[1]], "auto")
  expect_equal(project_status(f$con)$accepted_by, "auto")
})

test_that("missing depth parks instead of erroring", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = NA)
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
})

test_that("several failing rules are all reported", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3, topology = "linear")
  expect_setequal(codes(f), c("linear", "low_coverage"))
})

test_that("accept_linear passes a clean linear sample", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 80, topology = "linear")
  set_project_setting(f$con, "accept_linear", 1)
  expect_true(assembly_pass(f$con)$pass)
})

test_that("any curation warning parks; NA warnings parks", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE, warnings = c(0, 1, NA))
  expect_equal(annotation_pass(f$con)$pass, c(TRUE, FALSE, FALSE))
  expect_equal(annotation_pass(f$con)$reasons, c("", "warnings", "warnings"))
})

test_that("auto_accept never locks a parked unit; override does and logs the rules", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE, depth = 3)
  expect_equal(nrow(auto_accept(f$con, f$dir_out, "assemble")$accepted), 0)
  override_accept(f$con, f$dir_out, "assemble", data.frame(ID = "A"), reason = "re-sequencing not possible")
  lg <- DBI::dbGetQuery(f$con, "SELECT action, rules, reason FROM review_log")
  expect_equal(c(lg$action, lg$rules), c("override", "low_coverage"))
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_lock FROM assemble")[[1]], 1)
})

test_that("not_run and failed assembly codes", {
  f <- local_project_fixture(ids = c("A", "B"), assembled = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 3 WHERE ID = 'B'")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = NULL WHERE ID = 'A'")
  expect_equal(assembly_pass(f$con)$reasons, c("not_run", "failed"))
})

test_that("multiple paths and scaffolds count only non-ignored units", {
  u <- data.frame(ID = "A", path = c(1, 2, 2), scaffold = c(1, 1, 2))
  f <- local_project_fixture(ids = "A", annotate_units = u, assembled = TRUE)
  expect_setequal(codes(f), c("multiple_paths", "multiple_scaffolds"))
  DBI::dbExecute(f$con, "UPDATE assemblies SET ignore = 1 WHERE path = 2")
  expect_equal(codes(f), character())
})

test_that("each non-ignored unit needs enough depth", {
  u <- data.frame(ID = "A", path = 1, scaffold = c(1, 2))
  f <- local_project_fixture(ids = "A", annotate_units = u, assembled = TRUE, depth = c(80, 2))
  expect_equal(codes(f), c("multiple_scaffolds", "low_coverage"))
})

test_that("poor reference parks", {
  f <- local_project_fixture(ids = c("A", "B", "C"), assembled = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemble SET poor_blast_ref = 'poor' WHERE ID = 'A'")
  DBI::dbExecute(f$con, "UPDATE assemble SET poor_blast_ref = 'good' WHERE ID = 'B'")
  expect_equal(assembly_pass(f$con)$reasons, c("poor_reference", "", ""))
})

test_that("annotation not_run, failed and problematic codes", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 1 WHERE ID = 'A'")
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 3 WHERE ID = 'B'")
  DBI::dbExecute(f$con, "UPDATE annotate SET problematic = 'yes' WHERE ID = 'C'")
  expect_equal(annotation_pass(f$con)$reasons, c("not_run", "failed", "problematic"))
})

test_that("annotation_pass skips samples that are not assembly-locked", {
  f <- local_project_fixture(ids = c("A", "B"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 0 WHERE ID = 'B'")
  expect_equal(annotation_pass(f$con)$ID, "A")
})

test_that("auto_accept parks a stale sample as output_missing", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = "A", assembled = TRUE)
  r <- auto_accept(f$con, f$dir_out, "assemble")
  expect_equal(r$accepted$ID, "A")
  expect_equal(r$parked$reasons[r$parked$ID == "B"], "output_missing")
})

test_that("auto_accept on annotate parks a multi-path sample but locks the rest", {
  u <- data.frame(ID = c("A", "A", "B"), path = c(1, 2, 1), scaffold = 1)
  f <- local_project_fixture(ids = c("A", "B"), annotate_units = u, annotated = TRUE)
  r <- auto_accept(f$con, f$dir_out, "annotate")
  expect_equal(r$accepted$ID, "B")
  expect_equal(unique(r$parked$reasons[r$parked$ID == "A"]), "multiple_paths")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT ID FROM annotate WHERE annotate_lock = 1")$ID, "B")
})

test_that("accept_linear marks passing linear samples partial", {
  f <- local_project_fixture(ids = "A", on_disk = "A", annotated = TRUE, topology = "linear")
  set_project_setting(f$con, "accept_linear", 1)
  auto_accept(f$con, f$dir_out, "annotate")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT partial FROM annotate")[[1]], "yes")
})

test_that("override on annotate logs the failing codes", {
  f <- local_project_fixture(ids = "A", annotated = TRUE, warnings = 2)
  override_accept(f$con, f$dir_out, "annotate", data.frame(ID = "A", path = 1, scaffold = 1))
  expect_equal(DBI::dbGetQuery(f$con, "SELECT rules FROM review_log")[[1]], "warnings")
})

test_that("override reports lock blocks", {
  f <- local_project_fixture(ids = "A", assembled = TRUE)
  r <- override_accept(f$con, f$dir_out, "assemble", data.frame(ID = "A"))
  expect_equal(r$blocked$ID, "A")
  expect_equal(nrow(DBI::dbGetQuery(f$con, "SELECT * FROM review_log")), 0)
})

test_that("project_status reports stage, reasons, and manual locks", {
  f <- local_project_fixture(ids = c("A", "B", "C"), annotated = TRUE, warnings = c(0, 1, 0))
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 0, assemble_switch = 2 WHERE ID = 'C'")
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_lock = 1 WHERE ID = 'A'")
  assign_export_group(f$con, data.frame(ID = "A", path = 1, scaffold = 1), "g1")
  s <- project_status(f$con)
  expect_equal(s$stage, c("export", "annotate", "assemble"))
  expect_equal(s$accepted_by, c("manual", NA, NA))
  expect_equal(s$reasons, c("", "warnings", ""))
  expect_equal(s$export_group, c("g1", NA, NA))
})

test_that("project_status gives a sample with no annotate rows one row", {
  f <- local_project_fixture(ids = "A", assembled = TRUE)
  s <- project_status(f$con)
  expect_equal(nrow(s), 1)
  expect_true(is.na(s$path))
})

test_that("parse_stored_depth handles missing and blank", {
  expect_equal(parse_stored_depth("1 2 3"), 1:3)
  expect_length(parse_stored_depth(NA), 0)
  expect_length(parse_stored_depth(""), 0)
})

test_that("warning_tally counts semicolon separated warnings", {
  f <- local_project_fixture(ids = "A")
  DBI::dbAppendTable(f$con, "annotations", data.frame(
    ID = "A", path = 1, scaffold = 1, gene = c("g1", "g2"), pos1 = c(1, 2),
    warnings = c("x; y", "x")))
  expect_equal(warning_tally(f$con), data.frame(warning = c("x", "y"), n = c(2L, 1L)))
})

test_that("export_flags maps flag kinds to why codes and writes the report", {
  f <- local_project_fixture(ids = "A")
  proj <- list(con = f$con, db = f$con@dbname, dir_out = f$dir_out)
  flags <- dplyr::tibble(
    ID = c("A", "B"), gene = c("g1", "g2"),
    start_flag = c(TRUE, FALSE), stop_flag = c(FALSE, TRUE),
    identity_flag = c(TRUE, FALSE), internal_stop_flag = c(FALSE, FALSE))
  local_mocked_bindings(
    make_PCG_alignments = function(...) invisible(NULL),
    flag_PCG_outliers = function(...) list(
      flags = flags,
      internal_stops = data.frame(ID = "C", gene = "g3", n_stops = 1L)))
  r <- export_flags(proj, "g1")
  expect_equal(r$flagged$why, c("low_identity", "start_offset", "stop_offset", "internal_stop"))
  expect_equal(r$flagged$ID, c("A", "A", "B", "C"))
  expect_equal(basename(r$report), "AA_alignments_g1.html")
  expect_equal(dirname(r$report), file.path(f$dir_out, "export", "g1"))
})

test_that("switch codes outside 2 and 3 are not_run (allow-list)", {
  f <- local_project_fixture(ids = c("A", "B"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = 4 WHERE ID = 'A'")
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_switch = 4 WHERE ID = 'A'")
  expect_equal(annotation_pass(f$con)$reasons, c("not_run", ""))
  expect_equal(assembly_pass(f$con)$reasons, c("not_run", ""))
})

test_that("NA annotate_switch is not_run", {
  f <- local_project_fixture(ids = "A", annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE annotate SET annotate_switch = NULL")
  expect_equal(annotation_pass(f$con)$reasons, "not_run")
})

test_that("auto_accept locks nothing when dir_out is missing", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE)
  r <- auto_accept(f$con, file.path(f$path, "nope"), "assemble")
  expect_equal(nrow(r$accepted), 0)
  expect_equal(r$parked$reasons, "output_missing")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_lock FROM assemble")[[1]], 0)
})

test_that("override on unknown or unaccepted units is blocked and not logged", {
  f <- local_project_fixture(ids = c("A", "B"), annotated = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemble SET assemble_lock = 0 WHERE ID = 'B'")
  r <- override_accept(f$con, f$dir_out, "assemble", data.frame(ID = "ZZ"))
  expect_equal(r$blocked$reasons, "unknown_unit")
  u <- data.frame(ID = c("B", "Q"), path = 1, scaffold = 1)
  r <- override_accept(f$con, f$dir_out, "annotate", u)
  expect_setequal(r$blocked$reasons, c("assembly_not_accepted", "unknown_unit"))
  expect_equal(nrow(r$accepted), 0)
  expect_equal(nrow(DBI::dbGetQuery(f$con, "SELECT * FROM review_log")), 0)
  expect_equal(DBI::dbGetQuery(f$con, "SELECT annotate_lock FROM annotate WHERE ID = 'B'")[[1]], 0)
})

test_that("blocked samples are not marked partial", {
  f <- local_project_fixture(ids = c("A", "B"), on_disk = "A", assembled = TRUE, topology = "linear")
  DBI::dbExecute(f$con, "INSERT INTO annotate (ID, path, scaffold, partial) VALUES ('A', 1, 1, 'no'), ('B', 1, 1, 'no')")
  set_project_setting(f$con, "accept_linear", 1)
  auto_accept(f$con, f$dir_out, "assemble")
  expect_equal(DBI::dbGetQuery(f$con, "SELECT ID, partial FROM annotate ORDER BY ID")$partial, c("yes", "no"))
})

test_that("bad min_mean_depth setting gives low_coverage, never NA", {
  f <- local_project_fixture(ids = "A", assembled = TRUE)
  DBI::dbExecute(f$con, "UPDATE project_settings SET value = 'abc' WHERE key = 'min_mean_depth'")
  expect_equal(suppressWarnings(assembly_pass(f$con))$reasons, "low_coverage")
  f2 <- local_project_fixture(ids = "A", assembled = TRUE)
  DBI::dbExecute(f2$con, "UPDATE assemblies SET depth = '80 NA 80'")
  expect_equal(assembly_pass(f2$con)$reasons, "low_coverage")
})

test_that("sample with no or all-ignored assemblies, blank or NA depth, parks", {
  f <- local_project_fixture(ids = "A", assembled = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemblies SET depth = ''")
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
  DBI::dbExecute(f$con, "UPDATE assemblies SET depth = 'NA'")
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
  DBI::dbExecute(f$con, "UPDATE assemblies SET ignore = 1")
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
  DBI::dbExecute(f$con, "DELETE FROM assemblies")
  expect_equal(assembly_pass(f$con)$reasons, "low_coverage")
})

test_that("failed poor reference value parks", {
  f <- local_project_fixture(ids = "A", assembled = TRUE)
  DBI::dbExecute(f$con, "UPDATE assemble SET poor_blast_ref = 'failed'")
  expect_equal(assembly_pass(f$con)$reasons, "poor_reference")
})

test_that("a failed log write undoes the lock", {
  f <- local_project_fixture(ids = "A", on_disk = "A", assembled = TRUE)
  DBI::dbExecute(f$con, "DROP TABLE review_log")
  expect_error(auto_accept(f$con, f$dir_out, "assemble"))
  expect_equal(DBI::dbGetQuery(f$con, "SELECT assemble_lock FROM assemble")[[1]], 0)
})

test_that("export_flags uses the export defaults when thresholds are NULL", {
  f <- local_project_fixture(ids = "A")
  DBI::dbExecute(f$con, "UPDATE export_opts SET start_aa = NULL, stop_aa = NULL, ident_pct = NULL")
  proj <- list(con = f$con, db = f$con@dbname, dir_out = f$dir_out)
  got <- NULL
  local_mocked_bindings(
    make_PCG_alignments = function(group, db, out_path, start_aa, stop_aa, ident_pct) got <<- c(start_aa, stop_aa, ident_pct),
    flag_PCG_outliers = function(...) list(
      flags = dplyr::tibble(ID = character(), gene = character(), start_flag = logical(), stop_flag = logical(),
                            identity_flag = logical(), internal_stop_flag = logical()),
      internal_stops = data.frame(ID = character(), gene = character())))
  export_flags(proj, "g1")
  expect_equal(got, c(10, 10, 60))
})

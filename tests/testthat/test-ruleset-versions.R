# Each ruleset's settings are fingerprinted next to its version. If this fails,
# bump `ruleset_version` in the named params_<clade>_mito() and then refresh the
# file with: MITOPILOT_UPDATE_RULESET_FINGERPRINTS=1 devtools::test(filter = "ruleset-versions")
test_that("ruleset settings do not change without a version bump", {
  fp <- vapply(names(RULESET_MAP), function(t) {
    p <- do.call(paste0("params_", t), list())
    v <- p$ruleset_version
    p$ruleset_version <- NULL
    f <- withr::local_tempfile()
    writeLines(as.character(jsonlite::toJSON(p, auto_unbox = TRUE, digits = NA)), f)
    paste(v, unname(tools::md5sum(f)))
  }, character(1))
  path <- test_path("ruleset_fingerprints.csv")
  if (nzchar(Sys.getenv("MITOPILOT_UPDATE_RULESET_FINGERPRINTS"))) {
    utils::write.csv(data.frame(target = names(fp), fingerprint = unname(fp)),
                     path, row.names = FALSE)
  }
  old <- utils::read.csv(path)
  old <- stats::setNames(old$fingerprint, old$target)
  expect_setequal(names(fp), names(old))
  for (t in intersect(names(fp), names(old))) {
    same_version <- sub(" .*", "", fp[[t]]) == sub(" .*", "", old[[t]])
    if (same_version) {
      expect_identical(fp[[t]], old[[t]],
                       label = paste0(t, " settings (changed without a ruleset_version bump)"))
    } else {
      expect_identical(fp[[t]], old[[t]],
                       label = paste0(t, " fingerprint (version bumped; refresh ruleset_fingerprints.csv)"))
    }
  }
})

test_that("project ruleset label reports version, modification, and updates", {
  p <- params_fish_mito()
  expect_equal(ruleset_version_label("fish_mito", jsonlite::toJSON(p)), "v1.0")
  expect_equal(ruleset_version_label("fish_mito", jsonlite::toJSON(p, auto_unbox = TRUE)), "v1.0")
  p$ruleset_version <- NULL # project created before versioning
  expect_equal(ruleset_version_label("fish_mito", jsonlite::toJSON(p)), "v1.0")
  p$hit_threshold <- 1
  expect_equal(ruleset_version_label("fish_mito", jsonlite::toJSON(p)), "v1.0 (modified)")
  p$ruleset_version <- "0.9"
  expect_match(ruleset_version_label("fish_mito", jsonlite::toJSON(p)), "v0.9 \\(v1.0 available")
})

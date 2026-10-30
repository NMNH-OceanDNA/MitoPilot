# Extracted from test-synteny-js.R:44

# prequel ----------------------------------------------------------------------
.sy_state <- new.env()
sy_page <- function() {
  skip_if_not_installed("chromote")
  skip_if(!nzchar(Sys.getenv("CHROMOTE_CHROME")) && is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "no Chrome for chromote")
  if (isTRUE(.sy_state$chrome_failed)) skip("Chrome did not start")
  b <- tryCatch(chromote::ChromoteSession$new(), error = function(e) e)
  if (inherits(b, "error")) {
    .sy_state$chrome_failed <- TRUE
    skip(paste("Chrome did not start:", conditionMessage(b)))
  }
  withr::defer(b$close(), envir = parent.frame())
  file_url <- function(p) {
    p <- normalizePath(p, winslash = "/", mustWork = TRUE)
    paste0(if (startsWith(p, "/")) "file://" else "file:///", p)
  }
  html <- readLines(testthat::test_path("synteny", "index.html"), warn = FALSE)
  for (f in c("seqviewer.js", "synteny.js")) {
    js_file <- system.file("app", "www", f, package = "MitoPilot")
    html <- sub(sprintf('src="[^"]*%s"', f), sprintf('src="%s"', file_url(js_file)), html)
  }
  page <- testthat::test_path("synteny", "index.generated.html")
  writeLines(html, page)
  withr::defer(unlink(page), envir = parent.frame())
  b$Page$navigate(file_url(page))
  Sys.sleep(1)
  if (!identical(js(b, "typeof window.mpsyn"), "object")) {
    stop("synteny.js did not load")
  }
  b
}
js <- function(b, code) b$Runtime$evaluate(code, returnByValue = TRUE)$result$value
L <- function(r) jsonlite::fromJSON(r, simplifyVector = FALSE)

# test -------------------------------------------------------------------------
b <- sy_page()
r <- js(b, "(function(){var g=window.mpsyn.geom; var s='AC-GT', rr='ACGT-';
    var sn=g.nongap(s), rn=g.nongap(rr);
    return JSON.stringify([sn, rn, g.sampleCol(3,sn), g.sampleCol(9,sn),
      g.refCol(1,rn,4,0,false), g.refCol(1,rn,4,1,false), g.refCol(1,rn,4,0,true),
      Array.from(g.classify(s,rr))]);})()")
expect_equal(L(r), list(list(1, 2, 4, 5), list(1, 2, 3, 4), 4, 5, 1, 4, 4, list(0, 0, 2, 2, 2)))

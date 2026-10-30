# Extracted from test-synteny-js.R:97

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
sy_payload <- function() paste0(
  "{id:'sy-canvas', input:'sy-pick', sampleLabel:'S1', refAcc:'NC_1', refLabel:'Fish', refRc:false,",
  " s:'-'.repeat(100) + 'ACGT'.repeat(50) + '-'.repeat(108), r:'ACGT'.repeat(102),",
  " refLength:408, rotation:0, version:1, selected:null,",
  " sampleFeatures:[{row:1,gene:'nad1',type:'PCG',pos1:10,pos2:120,dir:'+'},",
  "                 {row:2,gene:'trnF',type:'tRNA',pos1:190,pos2:30,dir:'-'}],",
  " refFeatures:[{gene:'nad1',type:'PCG',pos1:14,pos2:124,dir:'+'}]}")

# test -------------------------------------------------------------------------
b <- sy_page()
r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); var s0 = window.mpsyn.state('sy-canvas');
    window.__handlers.mpsyn_select({id:'sy-canvas', row:1}); var s1 = window.mpsyn.state('sy-canvas');
    window.mpsyn.whole('sy-canvas'); var s2 = window.mpsyn.state('sy-canvas');
    window.mpsyn.zoom('sy-canvas', 2); var s3 = window.mpsyn.state('sy-canvas');
    window.mpsyn.zoom('sy-canvas', 100); var s4 = window.mpsyn.state('sy-canvas');
    var h = window.mpsyn.hitTest('sy-canvas', 100, 30);
    window.__handlers.mpsyn(%s); var s5 = window.mpsyn.state('sy-canvas');
    return JSON.stringify({s0:s0, s1:s1, s2:s2, s3:s3, s4:s4, s5:s5, h:h ? h.f.gene : null});})()",
    sy_payload(), sy_payload()))
x <- L(r)
expect_equal(x$s0$L, 408)
expect_equal(x$s0$sampleLanes, 1)
expect_equal(x$s0$refLanes, 1)

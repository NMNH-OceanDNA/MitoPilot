# tests/testthat/test-synteny-js.R
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

test_that("column maps follow gaps, rotation and flip", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var s='AC-GT', rr='ACGT-';
    var sn=g.nongap(s), rn=g.nongap(rr);
    return JSON.stringify([sn, rn, g.sampleCol(3,sn), g.sampleCol(9,sn),
      g.refCol(1,rn,4,0,false), g.refCol(1,rn,4,1,false), g.refCol(1,rn,4,0,true),
      Array.from(g.classify(s,rr))]);})()")
  expect_equal(L(r), list(list(1, 2, 4, 5), list(1, 2, 3, 4), 4, 5, 1, 4, 4, list(0, 0, 2, 1, 2)))
})

test_that("features split at the ends and flip with the reference", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var rn=[1,2,3,4,5,6,7,8,9,10];
    var fwd=g.refFeature({pos1:2,pos2:4,dir:'+'},{rNongap:rn,refLength:10,rotation:0,rc:false});
    var rc=g.refFeature({pos1:2,pos2:4,dir:'+'},{rNongap:rn,refLength:10,rotation:0,rc:true});
    var straddle=g.refFeature({pos1:9,pos2:2,dir:'-'},{rNongap:rn,refLength:10,rotation:0,rc:true});
    return JSON.stringify([fwd, rc, straddle, g.pieces(straddle.c1,straddle.c2,10), g.pieces(3,5,10)]);})()")
  expect_equal(L(r), list(list(c1 = 2, c2 = 4, dir = "+"), list(c1 = 7, c2 = 9, dir = "-"),
                          list(c1 = 9, c2 = 2, dir = "+"),
                          list(list(9, 10), list(1, 2)), list(list(3, 5))))
})

test_that("positions under a column undo flip and rotation, null on a gap", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var s='AC-GT', rr='ACGT-';
    var sc=g.counts(s), rc=g.counts(rr);
    return JSON.stringify([g.sampleAt(3,s,sc), g.sampleAt(4,s,sc), g.refAt(5,rr,rc,4,0,false),
      g.refAt(1,rr,rc,4,2,false), g.refAt(1,rr,rc,4,0,true), g.refAt(2,rr,rc,4,1,true)]);})()")
  expect_equal(L(r), list(NULL, 3, NULL, 3, 4, 4))
})

test_that("identity is the match fraction over columns with both bases", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom; var cls=g.classify('ACGTA-','ACCTAA');
    var m=g.prefix(cls,0), both=g.prefixNot(cls,2);
    return JSON.stringify([g.identity(1,6,m,both), g.identity(1,2,m,both), g.identity(6,6,m,both)]);})()")
  expect_equal(L(r), list(0.8, 1, NULL))
})

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

test_that("identity window widens a narrow pixel bin to at least 10 columns, clamped", {
  b <- sy_page()
  r <- js(b, "(function(){var g=window.mpsyn.geom;
    return JSON.stringify([g.window(5,5,100), g.window(50,51,100), g.window(1,20,100), g.window(99,100,100)]);})()")
  expect_equal(L(r), list(list(1, 10), list(46, 55), list(1, 20), list(95, 100)))
})

# 408-column alignment: sample bases occupy columns 101-300, reference all 408.
sy_payload <- function() paste0(
  "{id:'sy-canvas', input:'sy-pick', sampleLabel:'S1', refAcc:'NC_1', refLabel:'Fish', refRc:false,",
  " s:'-'.repeat(100) + 'ACGT'.repeat(50) + '-'.repeat(108), r:'ACGT'.repeat(102),",
  " refLength:408, rotation:0, version:1, selected:null,",
  " sampleFeatures:[{row:1,gene:'nad1',type:'PCG',pos1:10,pos2:120,dir:'+'},",
  "                 {row:2,gene:'trnF',type:'tRNA',pos1:190,pos2:30,dir:'-'}],",
  " refFeatures:[{gene:'nad1',type:'PCG',pos1:14,pos2:124,dir:'+'}]}")

test_that("the viewer opens on the sample span, fits, zooms and keeps its view on a same-version reload", {
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
  expect_equal(x$s0$L, 408); expect_equal(x$s0$sampleLanes, 2); expect_equal(x$s0$refLanes, 1)
  expect_true(x$s0$viewStart <= 101 && x$s0$viewStart > 50); expect_gt(x$s0$ppc, x$s2$ppc)
  expect_equal(x$s1$selected, 1)
  expect_true(x$s1$viewStart <= 110 && x$s1$viewStart + 940 / x$s1$ppc >= 221)
  expect_equal(x$s2$viewStart, 1); expect_equal(x$s3$ppc, x$s2$ppc * 2); expect_equal(x$s4$ppc, 14)
  expect_equal(x$h, "nad1")
  expect_equal(x$s5$ppc, 14); expect_equal(x$s5$viewStart, x$s4$viewStart)
})

test_that("clicking a sample arrow selects its row; a reference arrow only fits", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); window.mpsyn.whole('sy-canvas');
    var v = window.mpsyn.state('sy-canvas');
    var xs = 60 + (150 - v.viewStart) * v.ppc, xr = 60 + (60 - v.viewStart) * v.ppc;
    window.mpsyn.click('sy-canvas', xs, 22 + 4 + 8); var a = window.mpsyn.state('sy-canvas').selected;
    var n0 = window.__inputs.length;
    window.mpsyn.click('sy-canvas', xr, v.height - 22 - 4 - 8); var s = window.mpsyn.state('sy-canvas');
    return JSON.stringify({a:a, n0:n0, n1:window.__inputs.length, name:window.__inputs[0].name,
      row:window.__inputs[0].value.row, ppcWhole:v.ppc, ppcAfter:s.ppc, selAfter:s.selected});})()", sy_payload()))
  x <- L(r)
  expect_equal(x$a, 1); expect_equal(x$n0, 1); expect_equal(x$n1, 1)
  expect_equal(x$name, "sy-pick"); expect_equal(x$row, 1)
  expect_gt(x$ppcAfter, x$ppcWhole); expect_equal(x$selAfter, 1)
})

test_that("the tooltip names gaps and never prints NaN", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s);
    return JSON.stringify([window.mpsyn.tip('sy-canvas', 1), window.mpsyn.tip('sy-canvas', 101),
                           window.mpsyn.tip('sy-canvas', 102)]);})()", sy_payload()))
  x <- jsonlite::fromJSON(r)
  expect_match(x[1], "sample gap"); expect_match(x[1], "ref 1 A")
  expect_match(x[2], "sample 1 A"); expect_match(x[2], "ref 101 A"); expect_match(x[2], "match$")
  expect_match(x[3], "sample 2 C"); expect_match(x[3], "ref 102 C")
  expect_false(any(grepl("NaN", x)))
})

test_that("a replaced canvas gets a fresh viewer", {
  b <- sy_page()
  r <- js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); var old = document.getElementById('sy-canvas');
    var nw = document.createElement('canvas'); nw.id = 'sy-canvas'; old.parentNode.replaceChild(nw, old);
    window.__handlers.mpsyn(%s); return JSON.stringify(window.mpsyn.state('sy-canvas').L);})()",
    sy_payload(), sy_payload()))
  expect_equal(jsonlite::fromJSON(r), 408)
})

test_that("a canvas replaced after the payload arrived is redrawn from the cached payload", {
  b <- sy_page()
  js(b, sprintf("(function(){ window.__handlers.mpsyn(%s); var old = document.getElementById('sy-canvas');
    var nw = document.createElement('canvas'); nw.id = 'sy-canvas'; old.parentNode.replaceChild(nw, old); })()", sy_payload()))
  Sys.sleep(0.5)
  r <- js(b, "JSON.stringify(window.mpsyn.state('sy-canvas'))")
  x <- L(r)
  expect_true(x$live); expect_equal(x$L, 408); expect_gt(x$height, 0)
  expect_equal(js(b, "parseInt(document.getElementById('sy-canvas').style.height, 10) > 0"), TRUE)
})

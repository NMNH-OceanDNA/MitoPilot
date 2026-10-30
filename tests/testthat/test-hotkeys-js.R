# tests/testthat/test-hotkeys-js.R
# Browser tests for inst/app/www/hotkeys.js against a fake-DOM fixture.
.hk_state <- new.env()

hk_page <- function() {
  skip_if_not_installed("chromote")
  skip_if(!nzchar(Sys.getenv("CHROMOTE_CHROME")) && is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "no Chrome for chromote")
  if (isTRUE(.hk_state$chrome_failed)) skip("Chrome did not start")
  b <- tryCatch(chromote::ChromoteSession$new(), error = function(e) e)
  if (inherits(b, "error")) {
    .hk_state$chrome_failed <- TRUE
    skip(paste("Chrome did not start:", conditionMessage(b)))
  }
  withr::defer(b$close(), envir = parent.frame())
  file_url <- function(p) {
    p <- normalizePath(p, winslash = "/", mustWork = TRUE)
    paste0(if (startsWith(p, "/")) "file://" else "file:///", p)
  }
  js_file <- system.file("app", "www", "hotkeys.js", package = "MitoPilot")
  html <- readLines(testthat::test_path("hotkeys", "index.html"), warn = FALSE)
  html <- sub('src="[^"]*hotkeys\\.js"', sprintf('src="%s"', file_url(js_file)), html)
  page <- testthat::test_path("hotkeys", "index.generated.html")
  writeLines(html, page)
  withr::defer(unlink(page), envir = parent.frame())
  b$Page$navigate(file_url(page))
  Sys.sleep(1)
  if (!identical(js(b, "typeof window.mpHotkeys"), "object")) {
    stop("hotkeys.js did not load from ", js_file)
  }
  js(b, "localStorage.removeItem('mp_hotkeys')")
  b
}
js <- function(b, code) b$Runtime$evaluate(code, returnByValue = TRUE)$result$value
clicks <- function(b) {
  r <- js(b, "JSON.stringify(window.__clicks)")
  if (identical(r, "[]")) character(0) else jsonlite::fromJSON(r)
}
press <- function(b, key, shift = FALSE, target = "document") {
  js(b, sprintf(
    "(function(){var t=%s; t.dispatchEvent(new KeyboardEvent('keydown',{key:'%s',shiftKey:%s,bubbles:true,cancelable:true})); return true;})()",
    target, key, tolower(shift)))
  invisible(clicks(b))
}

test_that("hotkeys are off by default and the pill turns them on", {
  b <- hk_page()
  expect_false(js(b, "window.mpHotkeys.enabled()"))
  expect_true(js(b, "document.getElementById('mp-hotkeys-on') !== null"))
  press(b, "r")
  expect_equal(js(b, "window.__clicks.length"), 0L)
  js(b, "window.mpHotkeys.setEnabled(true)")
  expect_equal(js(b, "localStorage.getItem('mp_hotkeys')"), "1")
  js(b, "window.mpHotkeys.bind('page:Assemble','r','Refresh',function(){return window.mpHotkeys.click('refresh');})")
  expect_equal(press(b, "r"), "refresh")
})

test_that("keyOf spells keys as the spec does", {
  b <- hk_page()
  r <- js(b, "JSON.stringify(['l','L',',','<',' ','ArrowDown','?'].map(function(k){
    return window.mpHotkeys.keyOf(new KeyboardEvent('keydown',{key:k,shiftKey:(k==='L'||k==='<'||k==='?')}));})
    .concat([window.mpHotkeys.keyOf(new KeyboardEvent('keydown',{key:'ArrowRight',shiftKey:true}))]))")
  expect_equal(jsonlite::fromJSON(r), c("l", "L", ",", "<", "Space", "ArrowDown", "?", "Shift+ArrowRight"))
})

test_that("typing targets, modifier combos and SweetAlerts are ignored", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  js(b, "window.mpHotkeys.bind('page:Assemble','r','Refresh',function(){return window.mpHotkeys.click('refresh');})")
  press(b, "r", target = "document.getElementById('free-text')")
  expect_equal(js(b, "window.__clicks.length"), 0L)
  js(b, "document.dispatchEvent(new KeyboardEvent('keydown',{key:'r',ctrlKey:true,bubbles:true}))")
  expect_equal(js(b, "window.__clicks.length"), 0L)
  js(b, "var s=document.createElement('div'); s.className='swal2-container'; document.body.appendChild(s)")
  press(b, "r")
  expect_equal(js(b, "window.__clicks.length"), 0L)
})

test_that("context follows the open modal, then the Step", {
  b <- hk_page()
  expect_equal(js(b, "window.mpHotkeys.context()"), "page:Assemble")
  js(b, "document.querySelector('#mode input[value=Export]').checked = true")
  expect_equal(js(b, "window.mpHotkeys.context()"), "page:Export")
  js(b, "document.getElementById('shiny-modal').innerHTML = '<div class=\"modal in\"><button id=\"x-update_state\"></button></div>'")
  expect_equal(js(b, "window.mpHotkeys.context()"), "state")
  js(b, "document.getElementById('shiny-modal').innerHTML = '<div class=\"modal in\"><button id=\"annotate-annotations-close\"></button></div>'")
  expect_equal(js(b, "window.mpHotkeys.context()"), "details")
  js(b, "document.getElementById('shiny-modal').innerHTML = '<div class=\"modal in\"></div>'")
  expect_equal(js(b, "window.mpHotkeys.context()"), "modal")
})

test_that("? opens a cheat sheet for the current context and Escape closes it", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  js(b, "window.mpHotkeys.bind('page:Assemble','r','Refresh table',function(){return true;})")
  press(b, "?", shift = TRUE)
  expect_false(js(b, "document.getElementById('mp-hotkeys-help').hidden"))
  expect_match(js(b, "document.getElementById('mp-hotkeys-help').textContent"), "Refresh table")
  expect_match(js(b, "document.getElementById('mp-hotkeys-help').textContent"), "Mouse only")
  press(b, "Escape")
  expect_true(js(b, "document.getElementById('mp-hotkeys-help').hidden"))
})

# ---- Annotate details modal ----
D <- "annotate-annotations-"

test_that("details modal: lifecycle, flags and actions click their buttons", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  for (k in c("q", "v", "r", "p", "l", "d", "m", "u", "g", "o")) press(b, k)
  press(b, "Q", TRUE); press(b, "L", TRUE); press(b, "D", TRUE)
  press(b, "ArrowRight", TRUE); press(b, "ArrowLeft", TRUE); press(b, "Backspace")
  expect_equal(clicks(b), paste0(D, c("close", "ID_verified", "reviewed", "problematic", "partial", "delete", "merge",
    "unjoin", "assign_gene", "auto_assign_orfs", "lock", "lock_only", "restore", "nav_next", "nav_prev", "back_to_review")))
})

test_that("details modal: panel keys toggle, Shift opens and scrolls, alignment goes through its summary", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "t"); expect_false(js(b, sprintf("document.getElementById('%sannotation_table_details').open", D)))
  press(b, "s"); expect_true(js(b, sprintf("document.getElementById('%sseqview-section').open", D)))
  press(b, "s"); expect_false(js(b, sprintf("document.getElementById('%sseqview-section').open", D)))
  press(b, "S", TRUE); expect_true(js(b, sprintf("document.getElementById('%sseqview-section').open", D)))
  press(b, "b"); expect_true(js(b, sprintf("document.getElementById('%sblast_synteny_details').open", D)))
  press(b, "a"); expect_true(js(b, "window.__clicks.indexOf('align-summary') >= 0"))
  expect_true(js(b, sprintf("document.getElementById('%salignment_div').open", D)))
  press(b, "i"); expect_true(js(b, sprintf("document.getElementById('%snotes').closest('details').open", D)))
  press(b, "I", TRUE)
  expect_equal(js(b, "document.activeElement.id"), paste0(D, "notes"))
})

test_that("details modal: Escape cancels an open merge panel, Enter confirms it, otherwise Enter is inert", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "Enter"); expect_equal(js(b, "window.__clicks.length"), 0L)
  js(b, sprintf("document.getElementById('%smerge_select_div').style.display=''", D))
  press(b, "Enter"); press(b, "Escape")
  expect_equal(clicks(b), paste0(D, c("confirm_merge", "cancel_merge")))
})

test_that("details modal: copy keys hit the selected row's buttons", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "y"); press(b, "Y", TRUE)
  expect_equal(clicks(b), c("fas2", "faa2"))
})

test_that("details modal: zoom and viewer buttons target the hovered viewer, default sequence", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "+", TRUE); press(b, "F", TRUE); press(b, "W", TRUE); press(b, "C", TRUE)
  js(b, sprintf("document.getElementById('%ssynteny_canvas').dispatchEvent(new MouseEvent('mouseover',{bubbles:true}))", D))
  press(b, "-"); press(b, "F", TRUE)
  press(b, "G", TRUE)
  expect_equal(clicks(b), c(paste0("zoom:", D, "seqview-canvas:2"), "fit", "whole", "copy",
                            paste0("syn:", D, "synteny_canvas:0.5"), "synfit"))
  expect_equal(js(b, "document.activeElement.id"), paste0(D, "seqview-goto"))
})

test_that("details modal: bindings survive a re-rendered modal (Review Focus 2)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "v"); js(b, "window.detailsModal(); window.__clicks = []"); press(b, "v")
  expect_equal(clicks(b), paste0(D, "ID_verified"))
})

test_that("details modal: arrows walk the annotation rows without wrapping", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "ArrowDown"); press(b, "ArrowUp"); press(b, "Home"); press(b, "End")
  expect_equal(clicks(b), c("row3", "row1", "row1", "row3"))
})

test_that("details modal: Shift+Up/Down cycle join segments around the active one", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "ArrowDown", TRUE); press(b, "ArrowUp", TRUE)
  expect_equal(clicks(b), c("seg3", "seg1"))
})

test_that("boundary keys route to codon or nucleotide buttons by wrapper visibility (Review Focus 3)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, ","); expect_equal(js(b, "window.__inputs.length"), 0L)
  js(b, sprintf("document.getElementById('%sedit_mode_ctrls').style.display=''", D))
  press(b, ","); press(b, "."); press(b, "<", TRUE); press(b, ">", TRUE)
  js(b, sprintf("document.getElementById('%sstart_search_ctrl').style.display='none'", D))
  press(b, ","); press(b, ">", TRUE)
  r <- jsonlite::fromJSON(js(b, "JSON.stringify(window.__inputs.map(function(i){return i.name+'='+i.value;}))"))
  expect_equal(r, paste0(D, c("start-add=one", "start-minus=one", "stop-minus=one", "stop-add=one", "rrna-5-out=one", "stop-add=one")))
})

test_that("edit-mode letters click their controls", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  js(b, sprintf("document.getElementById('%sedit_mode_ctrls').style.display=''", D))
  press(b, "e"); press(b, "E", TRUE); press(b, "z"); press(b, "5"); press(b, "3"); press(b, "P", TRUE)
  press(b, "f"); press(b, "X", TRUE); press(b, "c")
  expect_equal(clicks(b), c(paste0(D, c("edit_mode", "save_edits", "discard_edits", "toggle_partial_start",
    "toggle_partial_stop", "polyA_stop")), "reduce_align", "local_blast"))
  expect_equal(jsonlite::fromJSON(js(b, "JSON.stringify(window.__inputs.slice(-1)[0])"))$name, paste0(D, "single_codon"))
})

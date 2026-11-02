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
press <- function(b, key, shift = FALSE, target = "document", alt = FALSE, code = "") {
  js(b, sprintf(
    "(function(){var t=%s; t.dispatchEvent(new KeyboardEvent('keydown',{key:'%s',code:'%s',shiftKey:%s,altKey:%s,bubbles:true,cancelable:true})); return true;})()",
    target, key, code, tolower(shift), tolower(alt)))
  invisible(clicks(b))
}

test_that("hotkeys are off by default and the pill turns them on", {
  b <- hk_page()
  expect_false(js(b, "window.mpHotkeys.enabled()"))
  expect_true(js(b, "document.getElementById('mp-hotkeys-on') !== null"))
  expect_true(js(b, "document.getElementById('mp-hotkeys-pill').parentElement === document.body"))
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
  expect_equal(js(b, "window.mpHotkeys.keyOf(new KeyboardEvent('keydown',{key:'{',shiftKey:true}))"), "{")
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
  js(b, "window.showModal('<button id=\"x-update_state\"></button>')")
  expect_equal(js(b, "window.mpHotkeys.context()"), "state")
  js(b, "window.showModal('<button id=\"annotate-annotate-annotations-close\"></button>')")
  expect_equal(js(b, "window.mpHotkeys.context()"), "details")
  js(b, "window.showModal('<p>Work dirs</p>')")
  expect_equal(js(b, "window.mpHotkeys.context()"), "modal")
  js(b, "window.closeModal()")
  expect_equal(js(b, "window.mpHotkeys.context()"), "page:Export")
})

test_that("? opens a cheat sheet for the current context and Escape closes it", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  js(b, "window.mpHotkeys.bind('page:Assemble','r','Refresh table',function(){return true;})")
  js(b, "window.mpHotkeys.bind('page:Assemble','R','Run reports',function(){return true;})")
  press(b, "?", shift = TRUE)
  expect_false(js(b, "document.getElementById('mp-hotkeys-help').hidden"))
  expect_match(js(b, "document.getElementById('mp-hotkeys-help').textContent"), "Refresh table")
  expect_match(js(b, "document.getElementById('mp-hotkeys-help').textContent"), "Shift+r", fixed = TRUE)
  expect_false(grepl("Shift+R", js(b, "document.getElementById('mp-hotkeys-help').textContent"), fixed = TRUE))
  expect_false(grepl("Mouse only", js(b, "document.getElementById('mp-hotkeys-help').textContent"), fixed = TRUE))
  press(b, "Escape")
  expect_true(js(b, "document.getElementById('mp-hotkeys-help').hidden"))
})

# ---- Annotate details modal ----
D <- "annotate-annotate-annotations-"

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

test_that("details modal: Shift+Up/Down walk the annotation rows without wrapping", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "ArrowDown", TRUE); press(b, "ArrowUp", TRUE); press(b, "Home"); press(b, "End")
  expect_equal(clicks(b), c("row3", "row1", "row1", "row3"))
})

test_that("details modal: plain Up/Down scroll the window and never touch rows", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "ArrowDown"); press(b, "ArrowDown")
  top <- js(b, "document.querySelector('#shiny-modal .modal-body').scrollTop")
  expect_gt(top, 0)
  press(b, "ArrowUp")
  expect_lt(js(b, "document.querySelector('#shiny-modal .modal-body').scrollTop"), top)
  expect_equal(clicks(b), character(0))
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
    "toggle_partial_stop", "polyA_stop")), "reduce_align"))
  expect_equal(jsonlite::fromJSON(js(b, "JSON.stringify(window.__inputs.slice(-1)[0])"))$name, paste0(D, "single_codon"))
})

# ---- sample tables and toolbar ----
test_that("table cursor skips hidden rows, Space selects, Enter opens details (Review Focus 4)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "ArrowDown"); press(b, "ArrowDown"); press(b, "ArrowDown")
  expect_equal(js(b, "document.querySelectorAll('#assemble-table [data-mp-kb-cursor]').length"), 1L)
  expect_equal(js(b, "document.querySelector('#assemble-table [data-mp-kb-cursor] input').getAttribute('onclick')"), "__clicks.push('sel4')")
  press(b, " "); press(b, "Enter"); press(b, "ArrowUp"); press(b, "x"); press(b, "O", TRUE); press(b, "ArrowUp"); press(b, "S", TRUE)
  expect_equal(clicks(b), c("sel4", "det4", "sel2", "out2", "spec1"))
  press(b, "Escape")
  expect_equal(js(b, "document.querySelectorAll('[data-mp-kb-cursor]').length"), 0L)
})

test_that("Shift+X extends the checked range from the last toggled row to the cursor", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "ArrowDown"); press(b, " ")
  js(b, "document.querySelector('#assemble-table .rt-tr input').checked = true")
  press(b, "ArrowDown"); press(b, "ArrowDown"); press(b, "X", TRUE)
  expect_equal(clicks(b), c("sel1", "sel2", "sel4"))
})

test_that("toolbar letters click the global buttons for the current Step", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "r"); press(b, "s"); press(b, "L", TRUE); press(b, "u"); press(b, "v"); press(b, "m"); press(b, "w"); press(b, "R", TRUE)
  press(b, "2"); press(b, "v"); press(b, "l"); press(b, "m")
  press(b, "3"); press(b, "e"); press(b, "g"); press(b, "G", TRUE); press(b, "a"); press(b, "]"); press(b, "[")
  expect_equal(clicks(b), c("refresh", "state", "lock", "run_modal", "assemble-meta_view_open", "workdir_browser-open_browser",
    "run_reports-open", "id_verified_top", "partial_top", "annotate-meta_view_open", "export", "group", "clear_group",
    "all:export-table", "page:3", "page:1"))
  expect_equal(js(b, "window.mpHotkeys.context()"), "page:Export")
})

test_that("/ focuses the active table's search box", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "/")
  expect_equal(js(b, "document.activeElement.className"), "rt-search")
})

test_that("state modal: digits pick a radio and Enter applies", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.stateModal()")
  press(b, "2"); press(b, "Enter")
  expect_equal(clicks(b), c("st2", "assemble-update_state"))
})

# ---- viewer panning ----
test_that("details modal: Left/Right pan the hovered viewer, even with its section closed (Review Focus 5)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  expect_false(js(b, sprintf("document.getElementById('%sseqview-section').open", D)))
  ok <- js(b, "(function(){ try { document.dispatchEvent(new KeyboardEvent('keydown',{key:'ArrowRight',bubbles:true,cancelable:true})); return true; } catch (e) { return String(e); } })()")
  expect_true(ok)
  press(b, "ArrowLeft")
  js(b, sprintf("document.getElementById('%ssynteny_canvas').dispatchEvent(new MouseEvent('mouseover',{bubbles:true}))", D))
  press(b, "ArrowRight")
  expect_equal(clicks(b), c(paste0("pan:", D, "seqview-canvas:0.25"), paste0("pan:", D, "seqview-canvas:-0.25"),
                            paste0("synpan:", D, "synteny_canvas:0.25")))
})

# ---- Assemble details and outlier review ----
test_that("assemble details: block nav, resolve mode, cursor and viewer keys", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.coverageModal()")
  press(b, "]"); press(b, "p"); press(b, "3"); press(b, "4"); press(b, "a"); press(b, "L", TRUE); press(b, "B", TRUE); press(b, "T", TRUE)
  press(b, "ArrowDown"); press(b, "ArrowDown"); press(b, "ArrowDown"); press(b, " "); press(b, "X", TRUE)
  press(b, "W", TRUE); press(b, "c"); press(b, "t"); press(b, "R", TRUE); press(b, "ArrowRight"); press(b, "G", TRUE)
  expect_equal(js(b, "document.activeElement.id"), "assemble-assemble-coverage_details-maptoref_viewer-goto")
  press(b, "q", target = "document.body")
  expect_equal(clicks(b), c(paste0("assemble-assemble-coverage_details-", c("next_block", "prev_block")), "res:iupac", "res:nmask",
    paste0("assemble-assemble-coverage_details-", c("align", "lock", "build_resolved", "trim_consensus")),
    "scaf2", "ignore2", "mtwhole", "cov", "nt", "reads",
    "pan:assemble-assemble-coverage_details-maptoref_viewer-canvas:0.25", "assemble-assemble-coverage_details-close_modal"))
})

test_that("assemble details: I focuses the notes box", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.coverageModal()")
  press(b, "I", TRUE)
  expect_equal(js(b, "document.activeElement.id"), "assemble-assemble-coverage_details-notes")
})

test_that("outlier review: gene arrows, sample cursor, per-row actions, no key for Continue", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.reviewModal()")
  press(b, "ArrowRight"); press(b, "ArrowLeft"); press(b, "]"); press(b, "[")
  press(b, "ArrowDown"); press(b, "h"); press(b, "x"); press(b, "Enter"); press(b, "S", TRUE); press(b, "Escape")
  expect_equal(clicks(b), c("export-review_next", "export-review_prev", "export-review_next", "export-review_prev",
                            "pick1", "resolved1", "goto1", "export-skip_gene"))
})


# ---- final review fixes ----
test_that("an ordinary modal blocks page keys (review Critical 1)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.showModal('<p>Work dirs</p>')")
  press(b, "r"); press(b, "L", TRUE); press(b, "ArrowDown"); press(b, " ")
  expect_equal(clicks(b), character(0))
  expect_equal(js(b, "document.querySelectorAll('[data-mp-kb-cursor]').length"), 0L)
})

test_that("duplicate toolbar ids: the visible copy is clicked (review Important 2)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "s"); press(b, "L", TRUE)
  expect_equal(clicks(b), c("state", "lock"))
})

test_that("the cursor marker survives a React className rewrite (review Important 4)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "ArrowDown")
  js(b, "document.querySelector('#assemble-table [data-mp-kb-cursor]').className = 'rt-tr rt-tr-selected'")
  press(b, " ")
  expect_equal(clicks(b), "sel1")
})

test_that("focus on a checkbox or the pill does not swallow keys; Escape leaves a text field (review Important 5)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "r", target = "document.getElementById('mp-hotkeys-on')")
  expect_equal(clicks(b), "refresh")
  press(b, "ArrowDown")
  press(b, " ", target = "document.querySelector('#assemble-table input.rt-select-input')")
  expect_equal(clicks(b), "refresh")
  js(b, "document.getElementById('free-text').focus()")
  press(b, "Escape", target = "document.getElementById('free-text')")
  expect_equal(js(b, "document.activeElement.tagName"), "BODY")
})

test_that("hovering the reads track still drives the main viewer (review Important 7)", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.coverageModal()")
  js(b, "document.getElementById('assemble-assemble-coverage_details-maptoref_viewer-reads').dispatchEvent(new MouseEvent('mouseover',{bubbles:true}))")
  press(b, "+", TRUE)
  expect_equal(clicks(b), "zoom:assemble-assemble-coverage_details-maptoref_viewer-canvas:2")
})


test_that("held arrow keys move the table cursor three rows per repeat", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true)")
  press(b, "ArrowDown")
  js(b, "document.dispatchEvent(new KeyboardEvent('keydown',{key:'ArrowDown',repeat:true,bubbles:true,cancelable:true}))")
  press(b, " ")
  expect_equal(clicks(b), "sel4")
})


test_that("zoom aliases: + or = in, - or _ out; the cheat sheet lists each pair once", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "="); press(b, "+", TRUE); press(b, "-"); press(b, "_", TRUE)
  expect_equal(clicks(b), paste0("zoom:", D, "seqview-canvas:", c(2, 2, 0.5, 0.5)))
  press(b, "?", TRUE)
  txt <- js(b, "document.getElementById('mp-hotkeys-help').textContent")
  expect_match(txt, "+ or =", fixed = TRUE)
  expect_match(txt, "- or _", fixed = TRUE)
  expect_equal(js(b, "Array.from(document.querySelectorAll('#mp-hotkeys-help kbd')).filter(function(k){return /[+=_-]/.test(k.textContent) && k.textContent.length <= 6;}).length"), 2L)
})

test_that("Left/Right scroll the alignment when it is hovered", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  press(b, "a")
  js(b, sprintf("document.querySelector('#%salignment_div .biojs_msa_rheader').dispatchEvent(new MouseEvent('mouseover',{bubbles:true}))", D))
  press(b, "ArrowRight"); press(b, "ArrowRight"); press(b, "ArrowLeft")
  expect_equal(js(b, sprintf("document.querySelector('#%salignment_div .biojs_msa_rheader').scrollLeft", D)), 25)
  press(b, "+", TRUE)
  expect_equal(clicks(b), "align-summary")
})


test_that("bracket keys move the boundary by the step box instead of one", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  js(b, sprintf("document.getElementById('%sedit_mode_ctrls').style.display=''", D))
  press(b, "["); press(b, "]"); press(b, "{", TRUE); press(b, "}", TRUE)
  r <- jsonlite::fromJSON(js(b, "JSON.stringify(window.__inputs.map(function(i){return i.name+'='+i.value;}))"))
  expect_equal(r, paste0(D, c("start-add=box", "start-minus=box", "stop-minus=box", "stop-add=box")))
  press(b, ",", alt = TRUE, code = "Comma"); press(b, "v", alt = TRUE, code = "KeyV")
  expect_equal(js(b, "window.__inputs.length"), 4L)
  expect_equal(clicks(b), character(0))
})



test_that("Shift+L unlocks when the window is locked (no lock button) in both details windows", {
  b <- hk_page(); js(b, "window.mpHotkeys.setEnabled(true); window.detailsModal()")
  js(b, sprintf("document.getElementById('%slock_only').remove()", D))
  press(b, "L", TRUE)
  expect_equal(clicks(b), paste0(D, "unlock"))
  js(b, "window.__clicks = []; window.coverageModal(); document.getElementById('assemble-assemble-coverage_details-lock').remove()")
  press(b, "L", TRUE)
  expect_equal(clicks(b), "assemble-assemble-coverage_details-unlock")
})

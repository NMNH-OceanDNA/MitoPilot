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

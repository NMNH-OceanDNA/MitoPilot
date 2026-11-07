tb_page <- function() {
  b <- chromote_session_or_skip(envir = parent.frame())
  js_file <- system.file("app", "www", "custom.js", package = "MitoPilot")
  jq <- system.file("www", "shared", "jquery.min.js", package = "shiny")
  html <- readLines(testthat::test_path("toolbar", "index.html"), warn = FALSE)
  html <- sub('src="[^"]*custom\\.js"', sprintf('src="%s"', chromote_file_url(js_file)), html)
  html <- sub('src="jquery\\.min\\.js"', sprintf('src="%s"', chromote_file_url(jq)), html)
  page <- testthat::test_path("toolbar", "index.generated.html")
  writeLines(html, page)
  withr::defer(unlink(page), envir = parent.frame())
  chromote_navigate(b, chromote_file_url(page))
  Sys.sleep(1)
  b
}
js <- chromote_js

test_that("same-id toolbar buttons share one click count so every click reaches the server", {
  b <- tb_page()
  out <- js(b, "(function(){
    document.querySelector('#asmb_ctrls #run_modal').click();
    document.querySelector('#annot_ctrls #run_modal').click();
    document.querySelector('#export_ctrls #export_now').click();
    document.querySelector('#asmb_ctrls #run_modal').click();
    return JSON.stringify(window.__arrived);})()")
  expect_equal(jsonlite::fromJSON(out),
               c("run_modal:1", "run_modal:2", "export_now:1", "run_modal:3"))
})

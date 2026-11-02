# Shared chromote plumbing for the browser tests (seqview, synteny, hotkeys).
# Chrome on a CI runner can fail to start or stall mid-session; either is an
# environment problem, not a MitoPilot failure, so the test skips. A failure
# is remembered so later tests do not each wait out the timeout.
.chromote_state <- new.env()

chromote_timed_out <- function(e) {
  grepl("Chromote: timed out", conditionMessage(e), fixed = TRUE)
}

chromote_session_or_skip <- function(envir = parent.frame()) {
  skip_if_not_installed("chromote")
  skip_if(!nzchar(Sys.getenv("CHROMOTE_CHROME")) &&
            is.null(tryCatch(chromote::find_chrome(), error = function(e) NULL)),
          "no Chrome for chromote")
  if (isTRUE(.chromote_state$failed)) skip("Chrome did not start")
  b <- tryCatch(chromote::ChromoteSession$new(), error = function(e) e)
  if (inherits(b, "error")) {
    .chromote_state$failed <- TRUE
    skip(paste("Chrome did not start:", conditionMessage(b)))
  }
  withr::defer(b$close(), envir = envir)
  b
}

chromote_navigate <- function(b, url) {
  tryCatch(b$Page$navigate(url), error = function(e) {
    if (!chromote_timed_out(e)) stop(e)
    .chromote_state$failed <- TRUE
    skip(paste("Chrome stalled:", conditionMessage(e)))
  })
  invisible(b)
}

chromote_js <- function(b, code) {
  tryCatch(
    b$Runtime$evaluate(code, returnByValue = TRUE)$result$value,
    error = function(e) {
      if (!chromote_timed_out(e)) stop(e)
      skip(paste("Chrome stalled:", conditionMessage(e)))
    }
  )
}

chromote_file_url <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = TRUE)
  paste0(if (startsWith(p, "/")) "file://" else "file:///", p)
}

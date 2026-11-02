# A stalled Chrome (Chromote timeout) skips the test instead of failing it.
stalled <- list(
  Page = list(navigate = function(url) stop("Chromote: timed out waiting for response to command Page.navigate")),
  Runtime = list(evaluate = function(...) stop("Chromote: timed out waiting for response to command Runtime.evaluate"))
)
broken <- list(Runtime = list(evaluate = function(...) stop("ReferenceError: x is not defined")))

test_that("a Chromote timeout on evaluate or navigate becomes a skip", {
  expect_condition(chromote_js(stalled, "1"), class = "skip")
  old <- .chromote_state$failed
  withr::defer(.chromote_state$failed <- old)
  expect_condition(chromote_navigate(stalled, "about:blank"), class = "skip")
  expect_true(isTRUE(.chromote_state$failed))
})

test_that("other browser errors still fail", {
  expect_error(chromote_js(broken, "x"), "ReferenceError")
})

rt_rv <- function(run_tools = 0L) {
  list(
    updating = data.frame(ID = "S1", annotate_opts = "default", curate_opts = "default"),
    annotate_opts = data.frame(annotate_opts = "default", cpus = 4, memory = 8,
                               ref_dir = "ref", ref_db = "Chordata", run_tools = run_tools),
    curate_opts = data.frame(curate_opts = "default", params = "{}")
  )
}

rt_modal <- function(user_annot, run_tools = 0L) {
  testthat::local_mocked_bindings(showModal = function(ui, ...) ui)
  s <- list(ns = shiny::NS("x"), userData = list(user_annot = user_annot))
  as.character(shiny::withReactiveDomain(s, annotate_opts_modal(rt_rv(run_tools), session = s)))
}

tools_box_tag <- function(html) regmatches(html, regexpr('<div[^>]*id="x-tools_box"[^>]*>', html))

test_that("run_tools checkbox only renders for user-annotation projects", {
  html_ua <- rt_modal(TRUE)
  html_std <- rt_modal(FALSE)
  expect_match(html_ua, "Run annotation tools")
  expect_match(html_ua, 'id="x-run_tools_box"')
  expect_match(html_ua, 'id="[^"]*tools_box"')
  expect_no_match(html_std, "Run annotation tools")
  expect_no_match(html_std, "run_tools_box")
})

test_that("Tools box is hidden only when user_annot and run_tools is off", {
  expect_match(tools_box_tag(rt_modal(TRUE, 0L)), "shinyjs-hide")
  expect_no_match(tools_box_tag(rt_modal(TRUE, 1L)), "shinyjs-hide")
  expect_no_match(tools_box_tag(rt_modal(FALSE, 0L)), "shinyjs-hide")
})

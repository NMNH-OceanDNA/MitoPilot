test_that("Hydra setup runs only on Hydra when qsub is missing", {
  skip_if(nzchar(Sys.which("qsub")), "qsub on PATH here")
  ran <- FALSE
  local_mocked_bindings(hydra_setup = function() { ran <<- TRUE; invisible(TRUE) })
  local_mocked_bindings(is_hydra_cluster = function() FALSE)
  expect_false(ensure_hydra_setup())
  expect_false(ran)
  local_mocked_bindings(is_hydra_cluster = function() TRUE)
  expect_message(ensure_hydra_setup(), "Hydra detected")
  expect_true(ran)
})

test_that("Hydra setup still runs when qsub is visible but nextflow is not", {
  ran <- FALSE
  local_mocked_bindings(hydra_setup = function() { ran <<- TRUE; invisible(TRUE) })
  local_mocked_bindings(is_hydra_cluster = function() TRUE)
  local_mocked_bindings(detect_qsub = function() TRUE, nf_on_path = function() TRUE)
  expect_false(ensure_hydra_setup())
  expect_false(ran)
  local_mocked_bindings(nf_on_path = function() FALSE)
  expect_message(ensure_hydra_setup(), "Hydra detected")
  expect_true(ran)
})

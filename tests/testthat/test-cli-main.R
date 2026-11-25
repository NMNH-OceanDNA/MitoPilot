test_that("cli_parse_args splits positionals, values, and flags", {
  p <- cli_parse_args(c("assemble", "--no-resume", "--project", "/x", "--min-depth=5", "--submit"),
                      flags = c("resume", "submit"))
  expect_equal(p$pos, "assemble")
  expect_identical(p$opts, list(resume = FALSE, project = "/x", min_depth = "5", submit = TRUE))
})
test_that("tokens after -- are positional", {
  p <- cli_parse_args(c("a", "--", "--b"))
  expect_equal(p$pos, c("a", "--b"))
})
test_that("value flag without a value fails with status 1", {
  e <- tryCatch(cli_parse_args("--project"), mitopilot_cli_error = identity)
  expect_s3_class(e, "mitopilot_cli_error"); expect_equal(e$status, 1L)
})
test_that("value flag followed by another flag fails", {
  expect_error(cli_parse_args(c("--project", "--json")), class = "mitopilot_cli_error")
})
test_that("unknown --no-x fails", {
  expect_error(cli_parse_args("--no-bogus", flags = "resume"), class = "mitopilot_cli_error")
})
test_that("cli_main returns 1 and writes nothing for unknown command", {
  expect_equal(suppressMessages(cli_main("nonsense")), 1L)
  expect_equal(suppressMessages(cli_main("--help")), 0L)
})
test_that("cli_main dispatches and reports errors", {
  fn <- function(args) if (length(args)) cli_fail("bad", 3L) else 7L
  testthat::local_mocked_bindings(get0 = function(x, ...) if (x == "cli_cmd_zzz") fn else NULL, .package = "base")
  expect_equal(cli_main("zzz"), 7L)
  expect_equal(suppressMessages(cli_main(c("zzz", "x"))), 3L)
  fn1 <- function(args) cli_fail("bad opt")
  testthat::local_mocked_bindings(get0 = function(x, ...) fn1, .package = "base")
  expect_message(cli_main("zzz"), "Run 'mitopilot zzz --help' for usage.", fixed = TRUE)
})
test_that("cli_check_opts rejects an unknown option", {
  expect_error(cli_check_opts(list(project = "/x", bogus = TRUE), c("project", "json")),
               class = "mitopilot_cli_error")
})
test_that("cli_ask re-asks until a valid choice", {
  con <- textConnection(c("9", "2")); withr::local_options(MitoPilot.cli_input = con)
  expect_equal(suppressMessages(cli_ask("Pick", c(a = "A", b = "B"))), "b")
})
test_that("cli_confirm fails without --yes when not interactive", {
  withr::local_options(MitoPilot.cli_interactive = FALSE)
  expect_error(cli_confirm("Rerun 3 samples?"), class = "mitopilot_cli_error")
  expect_true(cli_confirm("Rerun 3 samples?", yes = TRUE))
})
test_that("cli_confirm reads y/n when interactive", {
  withr::local_options(MitoPilot.cli_interactive = TRUE)
  con <- textConnection(c("n", "y")); withr::local_options(MitoPilot.cli_input = con)
  expect_false(suppressMessages(cli_confirm("Go?")))
  expect_true(suppressMessages(cli_confirm("Go?")))
})

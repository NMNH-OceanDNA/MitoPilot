#' @noRd
cli_is_interactive <- function() {
  getOption("MitoPilot.cli_interactive", isatty(stdin()))
}

#' @noRd
cli_readline <- function(prompt) {
  message(prompt, appendLF = FALSE)
  con <- getOption("MitoPilot.cli_input", file("stdin"))
  line <- readLines(con, n = 1L)
  if (length(line) == 0L) cli_fail("No input available")
  line
}

#' @noRd
cli_ask <- function(prompt, choices) {
  repeat {
    message(prompt)
    message(paste0("  ", seq_along(choices), ") ", unname(choices), collapse = "\n"))
    ans <- suppressWarnings(as.integer(trimws(cli_readline("> "))))
    if (length(ans) == 1L && !is.na(ans) && ans >= 1L && ans <= length(choices)) {
      return(names(choices)[ans])
    }
    message("Enter a number from 1 to ", length(choices))
  }
}

#' @noRd
cli_confirm <- function(prompt, yes = FALSE, detail = NULL) {
  if (yes) return(TRUE)
  if (!cli_is_interactive()) {
    cli_fail(paste0(prompt, " (needs --yes)"), subclass = "mitopilot_needs_confirm", detail = detail)
  }
  tolower(trimws(cli_readline(paste0(prompt, " [y/N] ")))) %in% c("y", "yes")
}

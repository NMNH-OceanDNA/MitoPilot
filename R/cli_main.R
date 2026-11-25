#' @noRd
cli_fail <- function(msg, status = 1L, subclass = NULL, detail = NULL) {
  stop(structure(
    class = c(subclass, "mitopilot_cli_error", "error", "condition"),
    list(message = msg, call = NULL, status = as.integer(status), detail = detail)
  ))
}

#' @noRd
cli_parse_args <- function(args, flags = character()) {
  pos <- character()
  opts <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--") {
      pos <- c(pos, args[-seq_len(i)])
      break
    }
    if (!startsWith(a, "--")) {
      pos <- c(pos, a)
    } else {
      key <- sub("^--", "", a)
      val <- NULL
      if (grepl("=", key, fixed = TRUE)) {
        val <- sub("^[^=]*=", "", key)
        key <- sub("=.*$", "", key)
      }
      key <- gsub("-", "_", key)
      flag_names <- gsub("-", "_", flags)
      if (!is.null(val)) {
        if (key %in% flag_names) {
          v <- tolower(val)
          if (v %in% c("yes", "true", "1")) val <- TRUE
          else if (v %in% c("no", "false", "0")) val <- FALSE
          else cli_fail(paste0("Option --", gsub("_", "-", key), " takes yes or no, got '", val, "'"))
        } else if (startsWith(key, "no_") && sub("^no_", "", key) %in% flag_names) {
          cli_fail(paste0("Option ", sub("=.*$", "", a), " takes no value"))
        }
        opts[[key]] <- val
      } else if (key %in% flag_names) {
        opts[[key]] <- TRUE
      } else if (startsWith(key, "no_")) {
        nm <- sub("^no_", "", key)
        if (!nm %in% flag_names) cli_fail(paste0("Unknown option: ", a))
        opts[[nm]] <- FALSE
      } else {
        if (i == length(args) || startsWith(args[i + 1L], "--")) cli_fail(paste0("Option ", a, " needs a value"))
        i <- i + 1L
        opts[[key]] <- args[i]
      }
    }
    i <- i + 1L
  }
  list(pos = pos, opts = opts)
}

#' @noRd
cli_check_opts <- function(opts, known) {
  bad <- setdiff(names(opts), gsub("-", "_", known))
  if (length(bad)) {
    cli_fail(paste0("Unknown option(s): ", paste0("--", gsub("_", "-", bad), collapse = ", ")))
  }
  invisible(TRUE)
}

cli_usage <- function() {
  paste(
    "Usage: mitopilot <command> [options]",
    "",
    "Commands:",
    "  init  status  run  accept  hold  requeue  rules",
    "  group  metadata  export  auto  cluster",
    "",
    "Run 'mitopilot <command> --help' for command options.",
    sep = "\n"
  )
}

#' @noRd
cli_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L || args[1] %in% c("--help", "-h", "help")) {
    message(cli_usage())
    return(0L)
  }
  fn <- get0(paste0("cli_cmd_", gsub("-", "_", args[1])),
    envir = asNamespace("MitoPilot"), mode = "function"
  )
  if (is.null(fn)) {
    message("Unknown command: ", args[1], "\n\n", cli_usage())
    return(1L)
  }
  tryCatch(
    as.integer(fn(args[-1])),
    mitopilot_cli_error = function(e) {
      msg <- conditionMessage(e)
      if (e$status == 1L) msg <- paste0(msg, "\nRun 'mitopilot ", args[1], " --help' for usage.")
      message(msg)
      e$status
    },
    error = function(e) {
      message("Error: ", conditionMessage(e))
      1L
    }
  )
}

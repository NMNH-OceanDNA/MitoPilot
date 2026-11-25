#' Install the mitopilot command
#'
#' Writes the `mitopilot` shell command, bound to a specific Rscript, into a
#' directory on your PATH.
#'
#' @param dir Directory to write the command into.
#' @param rscript Rscript executable the command should run.
#' @param overwrite Replace an existing command.
#' @return The path of the command, invisibly.
#' @export
install_cli <- function(dir = "~/.local/bin",
                        rscript = file.path(R.home("bin"), "Rscript"),
                        overwrite = TRUE) {
  dir <- path.expand(dir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  src <- system.file("cli", "mitopilot", package = "MitoPilot")
  if (!nzchar(src)) src <- app_sys("cli", "mitopilot")
  path <- file.path(dir, "mitopilot")
  if (file.exists(path) && !overwrite) {
    stop("File exists: ", path, call. = FALSE)
  }
  txt <- sub("@RSCRIPT@", rscript, readLines(src), fixed = TRUE)
  writeLines(txt, path)
  Sys.chmod(path, "0755")
  on_path <- normalizePath(dir) %in%
    normalizePath(strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]], mustWork = FALSE)
  message("Installed ", path)
  if (!on_path) {
    message("Add to PATH: export PATH=\"", dir, ":$PATH\"")
  }
  invisible(path)
}

local_test_db <- function(env = parent.frame()) {
  td <- withr::local_tempdir(.local_envir = env)
  mapping <- file.path(td, "mapping.csv")
  utils::write.csv(
    data.frame(ID = "S1", Taxon = "Danio rerio",
               R1 = "S1_R1.fastq.gz", R2 = "S1_R2.fastq.gz"),
    mapping, row.names = FALSE
  )
  db <- file.path(td, ".sqlite")
  new_db(db_path = db, mapping_fn = mapping,
         fetch_geome = FALSE, fetch_gbif = FALSE, fetch_ncbi = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  withr::defer(DBI::dbDisconnect(con), envir = env)
  con
}

local_project_fixture <- function(ids = "A", on_disk = character(),
                                  annotate_units = NULL, executor = "local",
                                  assembled = FALSE, depth = 80, topology = "circular",
                                  annotated = FALSE, warnings = 0,
                                  env = parent.frame()) {
  con <- local_test_db(env)
  path <- dirname(con@dbname)
  cfg <- sub("executor = 'local'", sprintf("executor = '%s'", executor),
             readLines(app_sys("config.local")), fixed = TRUE)
  writeLines(cfg, file.path(path, ".config"))
  dir_out <- file.path(path, "out")
  dir.create(dir_out)
  tmpl <- DBI::dbGetQuery(con, "SELECT * FROM assemble")
  ann <- DBI::dbGetQuery(con, "SELECT * FROM annotate")
  DBI::dbExecute(con, "DELETE FROM assemble")
  DBI::dbExecute(con, "DELETE FROM annotate")
  asm <- tmpl[rep(1, length(ids)), ]
  asm$ID <- ids
  if (annotated) assembled <- TRUE
  if (assembled) {
    asm$assemble_switch <- 2
    asm$topology <- rep_len(topology, length(ids))
  }
  if (annotated) asm$assemble_lock <- 1
  DBI::dbAppendTable(con, "assemble", asm)
  units <- annotate_units %||% data.frame(ID = ids, path = 1, scaffold = 1)
  asmbl <- data.frame(units, ignore = 0)
  if (assembled) {
    d <- rep_len(depth, nrow(units))
    asmbl$topology <- rep_len(topology, nrow(units))
    asmbl$depth <- ifelse(is.na(d), NA_character_, vapply(d, function(x) paste(rep(x, 5), collapse = " "), ""))
  }
  DBI::dbAppendTable(con, "assemblies", asmbl)
  if (annotated && is.null(annotate_units)) annotate_units <- units
  if (!is.null(annotate_units)) {
    a <- ann[rep(1, nrow(units)), ]
    a[c("ID", "path", "scaffold")] <- units[c("ID", "path", "scaffold")]
    if (annotated) {
      a$annotate_switch <- 2
      a$warnings <- rep_len(warnings, nrow(a))
      a$problematic <- "no"
      a$partial <- "no"
    }
    DBI::dbAppendTable(con, "annotate", a)
  }
  for (u in which(units$ID %in% on_disk)) {
    d <- file.path(dir_out, units$ID[u], "assemble", "default")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(d, paste0(units$ID[u], "_assembly_", units$path[u], ".fasta")))
  }
  list(con = con, path = path, dir_out = dir_out)
}

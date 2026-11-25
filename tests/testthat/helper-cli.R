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

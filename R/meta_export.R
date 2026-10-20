.meta_combos <- function(prefix) {
  switch(tolower(prefix), geome = GEOME_COMBOS, gbif = GBIF_COMBOS, ncbi = NCBI_COMBOS, NULL)
}

.meta_ymd <- function(y, m, d) {
  if (is.na(y) || !grepl("^[0-9]{4}$", y)) return(NA_character_)
  m <- suppressWarnings(as.integer(m))
  d <- suppressWarnings(as.integer(d))
  if (is.na(m) || m < 1 || m > 12) return(y)
  if (is.na(d) || d < 1 || d > 31) return(sprintf("%s-%02d", y, m))
  sprintf("%s-%02d-%02d", y, m, d)
}

.meta_key_col <- function(key) {
  p <- strsplit(key, ":", fixed = TRUE)[[1]]
  nm <- if (p[2] == "combo") p[3] else paste(p[3], paste(p[-(1:3)], collapse = ":"), sep = "_")
  paste0(tolower(p[1]), "_", gsub("[^A-Za-z0-9_]", "_", nm))
}

.meta_key_value <- function(recs, key) {
  p <- strsplit(key, ":", fixed = TRUE)[[1]]
  if (p[2] == "combo") {
    spec <- .meta_combos(p[1])[[p[3]]]
    return(if (is.null(spec)) NA_character_ else spec$fn(recs))
  }
  field <- paste(p[-(1:3)], collapse = ":")
  v <- recs$value[recs$level == p[3] & recs$field == field]
  if (length(v)) v[1] else NA_character_
}

meta_export_cols <- function(con, ids = NULL, keys = NULL) {
  .meta_ensure_tables(con)
  keys <- keys %||% DBI::dbGetQuery(con, "SELECT key FROM meta_export_fields")$key
  if (!length(keys)) return(NULL)
  recs <- DBI::dbGetQuery(con, "SELECT ID, source, level, depth, field, value FROM meta_records")
  if (is.null(ids)) ids <- unique(recs$ID)
  out <- data.frame(ID = ids)
  by_src <- split(recs, recs$source)
  by_src_id <- list()
  for (k in keys) {
    p <- strsplit(k, ":", fixed = TRUE)[[1]]
    src <- toupper(p[1])
    r <- by_src[[src]] %||% recs[0, ]
    out[[.meta_key_col(k)]] <- if (p[2] == "combo") {
      if (is.null(by_src_id[[src]])) by_src_id[[src]] <- split(r, r$ID)
      by_id <- by_src_id[[src]]
      vapply(ids, function(i) .meta_key_value(by_id[[i]] %||% r[0, ], k), character(1),
             USE.NAMES = FALSE)
    } else {
      field <- paste(p[-(1:3)], collapse = ":")
      h <- r[r$level == p[3] & r$field == field, ]
      h$value[match(ids, h$ID)]
    }
  }
  out
}

.meta_join <- function(dat, con) {
  g <- meta_export_cols(con, ids = unique(dat$ID))
  if (is.null(g)) return(dat)
  out <- dplyr::left_join(dat, g, by = "ID")
  cols <- setdiff(names(g), "ID")
  out[cols] <- lapply(out[cols], function(x) ifelse(is.na(x), "", x))
  out
}

meta_field_summary <- function(con, source) {
  .meta_ensure_tables(con)
  prefix <- tolower(source)
  recs <- DBI::dbGetQuery(con, "SELECT ID, level, depth, field, value FROM meta_records WHERE source = ?",
                          params = list(source))
  sel <- DBI::dbGetQuery(con, "SELECT key FROM meta_export_fields")$key
  defs <- .meta_combos(prefix)
  by_id <- split(recs, recs$ID)
  combos <- do.call(rbind, lapply(names(defs), function(nm) {
    vals <- vapply(by_id, defs[[nm]]$fn, character(1), USE.NAMES = FALSE)
    data.frame(key = paste0(prefix, ":combo:", nm), kind = "combo", level = NA_character_,
               field = paste(defs[[nm]]$sources, collapse = " + "),
               n_samples = sum(!is.na(vals)),
               example = if (any(!is.na(vals))) vals[!is.na(vals)][1] else NA_character_)
  }))
  raw <- if (nrow(recs)) {
    key <- paste(recs$level, recs$depth, recs$field, sep = "\r")
    first <- !duplicated(key)
    g <- recs[first, c("level", "depth", "field", "value")]
    n_ids <- vapply(split(recs$ID, key), function(x) length(unique(x)), integer(1))
    g$n_samples <- unname(n_ids[key[first]])
    g <- g[order(-g$depth, g$level, g$field), ]
    data.frame(
      key = paste0(prefix, ":raw:", g$level, ":", g$field), kind = "raw", level = g$level,
      field = g$field, n_samples = g$n_samples, example = g$value
    )
  }
  out <- rbind(combos, raw)
  out <- unique(out[, c("key", "kind", "level", "field", "n_samples", "example")])
  out$n_samples <- as.integer(out$n_samples)
  out$col <- vapply(out$key, .meta_key_col, character(1), USE.NAMES = FALSE)
  out$selected <- out$key %in% sel
  rownames(out) <- NULL
  out
}

.meta_save_fields <- function(con, keys) {
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM meta_export_fields")
    if (length(keys)) DBI::dbAppendTable(con, "meta_export_fields", data.frame(key = unique(keys)))
  })
  invisible(keys)
}

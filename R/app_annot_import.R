#' TRUE when a user annotation project still needs its gene names confirmed
#' @noRd
annot_import_blocks_run <- function(session) {
  isTRUE(session$userData$user_annot) && !user_annot_confirmed(session$userData$con)
}

#' Header notice for unconfirmed gene names, NULL once confirmed
#' @noRd
annot_import_notice <- function(con) {
  if (user_annot_confirmed(con)) return(NULL)
  n <- annot_import_counts(con)
  sprintf("Gene names not confirmed: %d unmapped of %d. Pipeline runs are blocked.",
          n$unmapped, n$all)
}

#' Status chip counts for the gene-name review
#' @noRd
annot_import_counts <- function(con) {
  x <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS a, SUM(COALESCE(manual, 0) = 1) AS e,
                               SUM(gene IS NULL AND COALESCE(manual, 0) != 1) AS u FROM gene_name_map")
  e <- as.integer(x$e %|NA|% 0L)
  u <- as.integer(x$u %|NA|% 0L)
  list(all = as.integer(x$a), unmapped = u, edited = e, auto = as.integer(x$a) - u - e)
}

#' Gene-name review rows, unmapped first, or in the order of `keys`
#' @noRd
annot_import_data <- function(con, keys = NULL) {
  gm <- DBI::dbGetQuery(con, "SELECT raw_name, feature_type, gene, manual, n_samples, n_features
                              FROM gene_name_map")
  ua <- DBI::dbGetQuery(con, "SELECT DISTINCT raw_name, feature_type, ID, source_file, note
                              FROM user_annotations ORDER BY ID")
  k <- paste(gm$raw_name, gm$feature_type, sep = "\r")
  ku <- factor(paste(ua$raw_name, ua$feature_type, sep = "\r"), levels = k)
  agg <- function(x, sep) {
    vapply(split(x, ku), function(v) paste(sort(unique(v[!is.na(v) & nzchar(v)])), collapse = sep), "")
  }
  ext <- tolower(tools::file_ext(ua$source_file))
  ftype <- dplyr::coalesce(c(gff3 = "GFF3", gff = "GFF3", bed = "BED", tbl = "TBL")[ext], toupper(ext))
  label <- ifelse(nzchar(gm$raw_name), gm$raw_name,
                  sprintf("(unnamed %s)", ifelse(nzchar(gm$feature_type), gm$feature_type, "BED feature")))
  stem <- sub("\\.\\.\\.$", "", gm$raw_name)
  full <- vapply(split(ua$note, ku), function(v) v[!is.na(v)][1] %|NA|% "", "")
  full <- ifelse(endsWith(gm$raw_name, "...") & startsWith(full, stem), full, label)
  status <- ifelse(gm$manual %in% 1L, "edited", ifelse(is.na(gm$gene), "unmapped", "auto"))
  cls <- vapply(gm$gene, function(g) if (is.na(g)) "none" else if (g == "ctrl") "ctrl" else mito_gene_type(g), "")
  d <- data.frame(label = label, raw_name = gm$raw_name, title = full, ftype = unname(agg(ftype, "/")), class = unname(cls),
                  ids = unname(agg(ua$ID, ", ")), n_samples = gm$n_samples, n_features = gm$n_features,
                  gene = gm$gene, status = status, feature_type = gm$feature_type)
  o <- if (is.null(keys)) {
    order(match(status, c("unmapped", "edited", "auto")), tolower(d$label), d$feature_type)
  } else {
    match(keys, k)
  }
  d <- d[o, , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' Imported gene-name review table
#' @noRd
annot_import_table <- function(con, ns, d = annot_import_data(con)) {
  tbl <- ns("tbl")
  js <- htmlwidgets::JS
  sel_filter <- function(id, choices) function(values, name) {
    tags$select(class = "form-control input-sm mp-gn-filter", `data-tbl` = tbl, `data-col` = id,
                `aria-label` = paste("Filter", name),
                onchange = sprintf("Reactable.setFilter('%s', '%s', this.value || undefined)", tbl, id),
                tags$option(value = "", "All"), lapply(choices, function(o) tags$option(value = o, o)))
  }
  exact <- js("function(rows, id, v) { return rows.filter(function(r) {
    return String(r.values[id]).split('/').indexOf(v) >= 0; }); }")
  opts <- c("(drop)", MITO_IMPORT_GENE_CHOICES)
  opts <- jsonlite::toJSON(unname(Map(c, opts, ifelse(opts == "ctrl", "ctrl (control region)", opts))))
  first <- if (any(d$status == "unmapped")) "unmapped" else ""
  reactable::reactable(
    d, compact = TRUE, pagination = FALSE, height = 420, searchable = TRUE, striped = TRUE,
    resizable = TRUE, onClick = "expand", class = "mp-gn-tbl",
    details = reactable::colDef(html = TRUE, details = js(sprintf(
      "function(ri) { return mpGN.detail(ri, '%s', '%s'); }", tbl, ns("open")))),
    columns = list(
      label = reactable::colDef(
        name = "Raw name", width = 220, filterable = TRUE, html = TRUE,
        cell = js("function(ci) { var v = mpGN.esc(ci.value); return ci.row.title === ci.value ? v :
          '<abbr class=\"mp-abbr\" title=\"' + mpGN.esc(ci.row.title) + '\">' + v + '</abbr>'; }")),
      raw_name = reactable::colDef(show = FALSE),
      title = reactable::colDef(show = FALSE),
      ftype = reactable::colDef(name = "File type", width = 90, filterable = TRUE, filterMethod = exact,
                                filterInput = sel_filter("ftype", c("GFF3", "BED", "TBL"))),
      class = reactable::colDef(name = "Class", width = 80, filterable = TRUE, filterMethod = exact,
                                filterInput = sel_filter("class", c("PCG", "tRNA", "rRNA", "ctrl", "none"))),
      ids = reactable::colDef(show = FALSE),
      n_samples = reactable::colDef(
        name = "Samples", width = 80, align = "center", html = TRUE,
        cell = js("function(ci) { return '<span title=\"' + mpGN.esc(ci.row.ids) + '\">' + ci.value + '</span>'; }")),
      n_features = reactable::colDef(name = "Features", width = 80, align = "center"),
      gene = reactable::colDef(name = "Gene", width = 180, filterable = TRUE, html = TRUE,
                               cell = js(sprintf("function(ci) { return mpGN.select(ci, '%s', %s); }",
                                                 ns("pick"), opts))),
      status = reactable::colDef(
        name = "Status", width = 115, html = TRUE, filterable = TRUE, filterMethod = exact,
        filterInput = sel_filter("status", c("unmapped", "edited", "auto")),
        cell = rt_pill(c(auto = "success", edited = "info", unmapped = "warning"))),
      feature_type = reactable::colDef(show = FALSE)
    )
  ) |> htmlwidgets::onRender(sprintf("function() { mpGN.init('%s', '%s'); }", tbl, first))
}

#' Detail panel for one imported name: every feature carrying it
#' @noRd
annot_import_detail <- function(s) {
  yes <- function(v) ifelse(v %in% 1L, "yes", "")
  cols <- list(Sample = s$ID, Type = s$feature_type, Features = s$n_features, Start = s$pos1,
               End = s$pos2, Strand = s$direction, `5' partial` = yes(s$partial_start),
               `3' partial` = yes(s$partial_stop), Codon = s$codon_start, Segment = s$segment_group,
               Pseudogene = yes(s$pseudo), Note = s$note, File = s$source_file)
  cols <- lapply(cols, function(x) dplyr::coalesce(as.character(x), ""))
  tags$table(
    class = "table table-condensed table-striped",
    tags$thead(tags$tr(lapply(names(cols), tags$th))),
    tags$tbody(lapply(seq_len(nrow(s)), function(i) tags$tr(lapply(cols, function(x) tags$td(x[i])))))
  )
}

#' Status chips above the gene-name review table
#' @noRd
annot_import_chips <- function(con, ns) {
  n <- annot_import_counts(con)
  tbl <- ns("tbl")
  chip <- function(v, lbl, k) {
    tags$button(type = "button", class = "btn btn-default btn-sm mp-gn-chip", `data-tbl` = tbl,
                `data-v` = v, onclick = sprintf("mpGN.chip('%s', '%s')", tbl, v), sprintf("%s (%d)", lbl, k))
  }
  div(class = "mp-gn-chips", role = "group", `aria-label` = "Filter by status",
      chip("", "All", n$all), chip("unmapped", "Unmapped", n$unmapped),
      chip("edited", "Edited", n$edited), chip("auto", "Auto", n$auto),
      tags$script(HTML(sprintf("mpGN.mark('%s');", tbl))))
}

#' Imported gene-name review module
#' @noRd
annot_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ver <- reactiveVal(0L)
    bump <- function() ver(isolate(ver()) + 1L)
    cver <- reactiveVal(0L)
    keys <- NULL

    output$tbl <- reactable::renderReactable({
      ver()
      d <- annot_import_data(session$userData$con)
      keys <<- paste(d$raw_name, d$feature_type, sep = "\r")
      annot_import_table(session$userData$con, ns, d)
    })
    output$chips <- renderUI({
      ver()
      cver()
      annot_import_chips(session$userData$con, ns)
    })

    show_modal <- function() {
      n <- DBI::dbGetQuery(session$userData$con, "
        SELECT (SELECT COUNT(*) FROM samples) AS s, (SELECT COUNT(*) FROM user_annotations) AS f")
      bump()
      modalDialog(
        title = mp_modal_title("Imported gene names", sprintf("%d samples, %d features", n$s, n$f)),
        size = "l",
        opts_help(
          "Pick the gene each imported name stands for, or (drop) to leave those features out. ",
          "Confirm writes the annotations and queues changed samples for annotation; ",
          "pipeline runs stay blocked until you confirm."
        ),
        uiOutput(ns("chips")),
        reactable::reactableOutput(ns("tbl")),
        tags$script(HTML(sprintf(
          "$('#shiny-modal').one('hidden.bs.modal', function() { Shiny.setInputValue('%s', Date.now(), {priority: 'event'}); });",
          ns("closed")))),
        footer = mp_footer(
          primary = actionButton(ns("confirm"), "Confirm"), dismiss = NULL,
          extra = actionButton(ns("reset"), "Reset guesses")
        )
      ) |> showModal()
    }

    init("annot_import")
    on("annot_import", show_modal())

    session$onFlushed(function() {
      if (!user_annot_confirmed(session$userData$con)) withReactiveDomain(session, isolate(show_modal()))
    }, once = TRUE)

    observeEvent(input$pick, {
      e <- input$pick
      set_gene_name_map(session$userData$con, e$raw_name, e$feature_type, e$gene)
      session$sendCustomMessage("mp-gn-data", list(
        tbl = ns("tbl"), data = lapply(annot_import_data(session$userData$con, keys), I)))
      cver(isolate(cver()) + 1L)
      trigger("annot_import_state")
    })
    observeEvent(input$open, {
      e <- input$open
      con <- session$userData$con
      s <- gene_name_map_samples(con, e$raw_name, e$feature_type)
      session$sendCustomMessage("mp-gn-detail", list(id = e$id, html = as.character(annot_import_detail(s))))
    })
    observeEvent(input$reset, {
      reset_gene_name_map(session$userData$con)
      bump()
      trigger("annot_import_state")
    })
    observeEvent(input$confirm, {
      res <- tryCatch(confirm_user_annotations(session$userData$con, session$userData$dir),
                      error = function(e) e)
      if (inherits(res, "error")) {
        mp_alert("Could not confirm gene names", conditionMessage(res), type = "error")
        return()
      }
      mp_toast(sprintf("%d features kept, %d dropped, %d sample(s) queued for annotation",
                       res$kept %|NA|% 0L, res$dropped %|NA|% 0L, length(res$requeued)),
               type = "success")
      for (w in res$warnings) mp_toast(w, type = "warning")
      removeModal()
      trigger("refresh_assemble")
      trigger("refresh_annotate")
      trigger("annot_import_state")
    })
    observeEvent(input$closed, {
      if (!user_annot_confirmed(session$userData$con)) {
        mp_toast("Gene names not confirmed; pipeline runs stay blocked", type = "warning")
      }
    })
  })
}

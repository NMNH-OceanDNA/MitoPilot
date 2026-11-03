#' TRUE when a user annotation project still needs its gene names confirmed
#' @noRd
annot_import_blocks_run <- function(session) {
  isTRUE(session$userData$user_annot) && !user_annot_confirmed(session$userData$con)
}

#' Header notice for unconfirmed gene names, NULL once confirmed
#' @noRd
annot_import_notice <- function(con) {
  if (user_annot_confirmed(con)) return(NULL)
  n <- DBI::dbGetQuery(con, "SELECT SUM(gene IS NULL) AS u, COUNT(*) AS t FROM gene_name_map")
  sprintf("Gene names not confirmed: %d unmapped of %d. Pipeline runs are blocked.",
          as.integer(n$u %|NA|% 0L), as.integer(n$t))
}

#' Imported gene-name review table
#' @noRd
annot_import_table <- function(con, ns) {
  gm <- DBI::dbGetQuery(con, "SELECT raw_name, feature_type, gene, manual, n_samples, n_features
                              FROM gene_name_map
                              ORDER BY gene IS NOT NULL, feature_type, raw_name")
  status <- ifelse(gm$manual %in% 1L, "edited", ifelse(is.na(gm$gene), "unmapped", "auto"))
  pick <- function(value, index) {
    cur <- gm$gene[index] %|NA|% "(drop)"
    opts <- unique(c("(drop)", MITO_IMPORT_GENE_CHOICES, cur))
    as.character(tags$select(
      class = "form-control input-sm", `aria-label` = paste("Gene for", gm$raw_name[index]),
      onchange = sprintf(
        "Shiny.setInputValue('%s', {raw_name: %s, feature_type: %s, gene: this.value}, {priority: 'event'})",
        ns("pick"), jsonlite::toJSON(gm$raw_name[index], auto_unbox = TRUE),
        jsonlite::toJSON(gm$feature_type[index], auto_unbox = TRUE)),
      lapply(opts, function(o) tags$option(value = o, selected = if (identical(o, cur)) NA,
                                           if (o == "ctrl") "ctrl (control region)" else o))
    ))
  }
  d <- data.frame(raw_name = gm$raw_name,
                  feature_type = ifelse(nzchar(gm$feature_type), gm$feature_type, "BED"),
                  n_samples = gm$n_samples, n_features = gm$n_features,
                  gene = gm$gene, status = status)
  reactable::reactable(
    d, compact = TRUE, pagination = FALSE,
    columns = list(
      raw_name = reactable::colDef(name = "Raw name", minWidth = 160),
      feature_type = reactable::colDef(name = "File type", width = 110),
      n_samples = reactable::colDef(name = "Samples", width = 90),
      n_features = reactable::colDef(name = "Features", width = 90),
      gene = reactable::colDef(name = "Gene", width = 140, html = TRUE, cell = pick),
      status = reactable::colDef(name = "Status", width = 100, html = TRUE,
                                 cell = rt_pill(c(auto = "success", edited = "info", unmapped = "warning")))
    )
  )
}

#' Imported gene-name review module
#' @noRd
annot_import_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ver <- reactiveVal(0L)
    bump <- function() ver(isolate(ver()) + 1L)

    output$tbl <- reactable::renderReactable({
      ver()
      annot_import_table(session$userData$con, ns)
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
      bump()
      trigger("annot_import_state")
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

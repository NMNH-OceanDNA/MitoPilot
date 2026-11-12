# Trial build of the Quarto site. Run from the package root: Rscript website/build.R
# Guides come from vignettes/, news from NEWS.md, and the function reference is
# rendered by pkgdown and then re-wrapped as Quarto pages.
w <- "website"

# stage sources
unlink(file.path(w, c("articles", "reference", "news", "man")), recursive = TRUE)
for (d in c("articles", "reference", "news", "man")) dir.create(file.path(w, d))
file.copy(list.files("vignettes", "\\.Rmd$", full.names = TRUE), file.path(w, "articles"))
file.copy("vignettes/figures", file.path(w, "articles"), recursive = TRUE)
# like pkgdown: a guide that uses level-1 section headings gets them demoted,
# so they sit below the page title
for (f in list.files(file.path(w, "articles"), "\\.Rmd$", full.names = TRUE)) {
  x <- readLines(f)
  in_code <- cumsum(grepl("^```", x)) %% 2 == 1
  if (any(grepl("^# ", x) & !in_code)) {
    x <- append(x, "shift-heading-level-by: 1", after = 1)
    writeLines(x, f)
  }
}
file.copy("man/figures", file.path(w, "man"), recursive = TRUE)
file.copy(c("pkgdown/extra.css", "pkgdown/favicon/favicon.ico",
            "pkgdown/assets/ruleset-browser.html"), w, overwrite = TRUE)
writeLines(c("---", "title: News", "---", "", "{{< include ../../NEWS.md >}}"),
           file.path(w, "news", "index.qmd"))

# function reference: pkgdown renders, we keep only the page content
tmp <- tempfile("pkgdown")
pkg <- pkgdown::as_pkgdown(".", override = list(destination = tmp))
pkgdown::init_site(pkg)
pkgdown::build_reference(pkg)

raw_html <- function(nodes) {
  if (length(nodes) == 0) return(character(0))
  c("````{=html}", vapply(nodes, as.character, ""), "````", "")
}
for (f in list.files(file.path(tmp, "reference"), "\\.html$", full.names = TRUE)) {
  main <- xml2::xml_find_first(xml2::read_html(f), "//main")
  title <- trimws(xml2::xml_text(xml2::xml_find_first(main, ".//h1")))
  xml2::xml_remove(xml2::xml_find_all(main, "./div[contains(@class, 'page-header')]"))
  body <- unlist(lapply(xml2::xml_children(main), function(n) {
    kids <- xml2::xml_children(n)
    if (length(kids) > 0 && xml2::xml_name(kids[[1]]) == "h2") {
      # real markdown headings so Quarto can build the page outline
      c(sprintf("## %s {#%s}", trimws(xml2::xml_text(kids[[1]])),
                xml2::xml_attr(kids[[1]], "id")), "", raw_html(kids[-1]))
    } else {
      raw_html(list(n))
    }
  }))
  writeLines(
    c("---", paste0("title: ", jsonlite::toJSON(title, auto_unbox = TRUE)), "---", "", body),
    file.path(w, "reference", sub("\\.html$", ".qmd", basename(f)))
  )
}

# sidebar grouped like the pkgdown reference index
idx <- xml2::xml_find_first(
  xml2::read_html(file.path(tmp, "reference", "index.html")), "//main")
sections <- list()
for (n in xml2::xml_children(idx)) {
  h <- xml2::xml_find_first(n, "./h2")
  if (!is.na(h)) {
    sections[[length(sections) + 1]] <- list(
      section = trimws(xml2::xml_text(h)), contents = list())
  }
  a <- xml2::xml_find_all(n, ".//dt//a[@href]")
  a <- a[!duplicated(xml2::xml_attr(a, "href"))]
  if (length(a) == 0 || length(sections) == 0) next
  i <- length(sections)
  sections[[i]]$contents <- c(sections[[i]]$contents, unname(Map(
    function(text, href) {
      list(text = text, href = file.path("reference", sub("\\.html$", ".qmd", href)))
    }, trimws(xml2::xml_text(a)), xml2::xml_attr(a, "href"))))
}
yaml::write_yaml(
  list(website = list(sidebar = list(list(
    id = "reference", title = "Reference", `collapse-level` = 1L,
    contents = c(list(list(text = "Index", href = "reference/index.qmd")), sections)
  )))),
  file.path(w, "reference", "_sidebar.yml")
)

stopifnot(system2("quarto", c("render", w)) == 0)

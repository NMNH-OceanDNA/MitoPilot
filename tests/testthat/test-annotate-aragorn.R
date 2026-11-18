test_that("ARAGORN batch output is parsed in plain and mitochondrial modes", {
  out_lines <- c(
    ">ctg1 some description",
    "3 genes found",
    "1   mtRNA-Phe                       [1,69]\t31  \t(gaa)",
    "2   mtRNA-Gln                 c[3895,3965]\t33  \t(ttg)",
    "3   tRNA-Leu                   [2771,2844]\t36  \t(taa)"
  )
  testthat::local_mocked_bindings(run_tool = function(cmd, args, ...) {
    writeLines(out_lines, args[which(args == "-o") + 1])
  })
  asm <- Biostrings::DNAStringSet(c(ctg1 = strrep("A", 5000)))
  res <- suppressMessages(annotate_aragorn(asm, genetic_code = 2))
  expect_equal(res$gene, c("trnF", "trnQ", "trnL"))
  expect_equal(res$pos1, c(1L, 3895L, 2771L))
  expect_equal(res$direction, c("+", "-", "+"))
  expect_equal(res$anticodon, c("GAA", "TTG", "TAA"))
})

test_that("ARAGORN defaults to the metazoan mitochondrial search", {
  expect_equal(formals(annotate_aragorn)$aragorn_opts, "-mtx")
})

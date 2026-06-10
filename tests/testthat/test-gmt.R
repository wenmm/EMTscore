test_that("write_gmt and read_gmt round-trip genes", {
  gene_sets <- list(
    E = c("CDH1", "EPCAM"),
    M = c("VIM", "ZEB1", "FN1")
  )
  f <- tempfile(fileext = ".gmt")
  write_gmt(gene_sets, f)
  expect_true(file.exists(f))

  genes <- read_gmt(f)
  expect_s3_class(genes, "data.frame")
  expect_true(all(c("CDH1", "VIM", "ZEB1") %in% genes$gene))
  expect_equal(nrow(genes), 5L)
})

test_that("write_gmt rejects unnamed lists", {
  expect_error(write_gmt(list(c("A", "B")), tempfile()))
})

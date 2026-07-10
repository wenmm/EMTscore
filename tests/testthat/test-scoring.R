# Tests for the single-gene-set scoring functions.
# Each Execute_* function should return a data.frame with one row per sample
# (120 for the bundled geneExp) and a single column named after `score_names`.

test_that("Execute_AUCell returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_AUCell(geneExp, gmt_file, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_GSVA returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_GSVA(geneExp, gmt_file, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_ssGSEA returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_ssGSEA(geneExp, gmt_file, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_SCSE returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_SCSE(geneExp, gmt_file, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_JAS returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_JAS(geneExp, gmt_file, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_nnPCA returns a 120-row data.frame with finite scores", {
  skip_if_not_installed("nsprcomp")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_nnPCA(geneExp, gmt_file, dimension = 1, score_names = "score")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "score")
  expect_true(all(is.finite(res$score)))
})

test_that("Execute_nnPCA errors when no genes match the expression matrix", {
  skip_if_not_installed("nsprcomp")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gene_sets <- list(NoMatch = c("NOT_A_GENE_1", "NOT_A_GENE_2", "NOT_A_GENE_3"))
  gmt_file <- tempfile(fileext = ".gmt")
  write_gmt(gene_sets, gmt_file)

  expect_error(Execute_nnPCA(geneExp, gmt_file, dimension = 1, score_names = "score"))
})

test_that("Execute_nnPCA supports multiple dimensions", {
  skip_if_not_installed("nsprcomp")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- Execute_nnPCA(geneExp, gmt_file,
    dimension = 2, score_names = c("M1_score", "M2_score")
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), c("M1_score", "M2_score"))
})

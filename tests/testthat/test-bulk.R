# Tests for add_EMT_score_Bulk (bulk EMT scoring wrapper).
# It accepts a named list of expression matrices and returns a named list of
# data.frames, each with one row per sample and a column named after `emt_name`.

test_that("add_EMT_score_Bulk (AUCell) returns a named list of score data.frames", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- add_EMT_score_Bulk(
    expr_mat_list = list(s = geneExp), gmt_file,
    emt_name = "EMT", method = "AUCell", dimension = 1
  )
  expect_type(res, "list")
  expect_named(res, "s")
  expect_s3_class(res$s, "data.frame")
  expect_equal(nrow(res$s), 120L)
  expect_true("EMT" %in% colnames(res$s))
  expect_true(all(is.finite(res$s$EMT)))
  # rownames should be the sample names from geneExp
  expect_setequal(rownames(res$s), colnames(geneExp))
})

test_that("add_EMT_score_Bulk (nnPCA) works on a single matrix", {
  skip_if_not_installed("nsprcomp")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- add_EMT_score_Bulk(
    expr_mat_list = list(s = geneExp), gmt_file,
    emt_name = "EMT", method = "nnPCA", dimension = 1
  )
  expect_s3_class(res$s, "data.frame")
  expect_equal(nrow(res$s), 120L)
  expect_true(all(is.finite(res$s$EMT)))
})

test_that("add_EMT_score_Bulk handles a named list of several matrices", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  res <- add_EMT_score_Bulk(
    expr_mat_list = list(A = geneExp, B = geneExp), gmt_file,
    emt_name = "EMT", method = "AUCell", dimension = 1
  )
  expect_named(res, c("A", "B"))
  expect_equal(nrow(res$A), 120L)
  expect_equal(nrow(res$B), 120L)
})

test_that("add_EMT_score_Bulk errors when gmt_file is NULL", {
  data(geneExp)
  expect_error(
    add_EMT_score_Bulk(
      expr_mat_list = list(s = geneExp), gmt_file = NULL,
      emt_name = "EMT", method = "AUCell", dimension = 1
    ),
    "gmt_file"
  )
})

test_that("add_EMT_score_Bulk rejects an unsupported method", {
  data(geneExp)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
  expect_error(
    add_EMT_score_Bulk(
      expr_mat_list = list(s = geneExp), gmt_file,
      emt_name = "EMT", method = "not_a_method", dimension = 1
    )
  )
})

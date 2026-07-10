# Tests for the parallel multi-gene-set scoring functions.
# EM_signature.gmt contains 2 gene sets. With cores = 1 these run quickly.
#
# Return-shape note (verified against the bundled data):
#   All *_parallel scorers (AUCell / GSVA / ssGSEA / JASMINE / nnPCA / SCSE)
#     return a MATRIX with samples as rows (120) and gene sets as columns (2).

gmt_multi <- function() {
  system.file("extdata", "EM_signature.gmt", package = "EMTscore")
}

expected_sets <- c("Panchy_et_al_E_signature", "Panchy_et_al_M_signature")

test_that("Execute_AUCell_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_AUCell_parallel(geneExp, gmt_multi(), cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

test_that("Execute_GSVA_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_GSVA_parallel(geneExp, gmt_multi(), cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

test_that("Execute_ssGSEA_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_ssGSEA_parallel(geneExp, gmt_multi(), cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

test_that("Execute_JASMINE_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_JASMINE_parallel(geneExp, gmt_multi(), cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

test_that("Execute_nnPCA_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("nsprcomp")
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_nnPCA_parallel(geneExp, gmt_multi(), dimension = 1, cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

test_that("Execute_SCSE_parallel returns a 120 x 2 matrix over gene sets", {
  skip_if_not_installed("GSA")
  skip_if_not_installed("BiocParallel")
  data(geneExp)

  res <- Execute_SCSE_parallel(geneExp, gmt_multi(), cores = 1)
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_setequal(colnames(res), expected_sets)
  expect_true(all(is.finite(res)))
})

# Tests for compute_Signature_score on a bulk expression matrix.
# For matrix input it returns a numeric matrix with one row per sample (120)
# and a single column named after `score_name`, min-max scaled to [0, 1].

test_that("compute_Signature_score returns a 120 x 1 scaled matrix", {
  data(geneExp)
  signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")

  res <- compute_Signature_score(geneExp, signature_file, "Stemness_Score")
  expect_true(is.matrix(res))
  expect_equal(nrow(res), 120L)
  expect_equal(colnames(res), "Stemness_Score")
  expect_true(all(is.finite(res)))
  # scores are min-max normalised to [0, 1]
  expect_equal(min(res), 0)
  expect_equal(max(res), 1)
  # rownames are the sample names
  expect_setequal(rownames(res), colnames(geneExp))
})

test_that("compute_Signature_score errors on a missing signature file", {
  data(geneExp)
  expect_error(
    compute_Signature_score(geneExp, tempfile(fileext = ".tsv"), "Score"),
    "does not exist"
  )
})

test_that("compute_Signature_score rejects non-matrix, non-Seurat input", {
  signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
  expect_error(
    compute_Signature_score(list(1, 2, 3), signature_file, "Score")
  )
})

# Tests for data_prepare, which merges a cell/sample annotation table with a
# score result (by row.names) and sets the annotation `merge_colname` as row
# names of the annotation prior to the merge.

test_that("data_prepare merges annotation with scores and sets rownames", {
  skip_if_not_installed("AUCell")
  skip_if_not_installed("GSEABase")
  data(geneExp)
  data(cell_annotation_file)
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  scores <- Execute_AUCell(geneExp, gmt_file, score_names = "EMT")
  merged <- data_prepare(cell_annotation_file, scores, merge_colname = "name")

  expect_s3_class(merged, "data.frame")
  # annotation columns and the score column are all present
  expect_true(all(c("name", "source", "celltype_annotation", "EMT") %in%
    colnames(merged)))
  # the score column carries finite numbers
  expect_true(all(is.finite(merged$EMT)))
  # rownames are set (inner join keeps only samples present in both inputs)
  expect_false(is.null(rownames(merged)))
  expect_true(nrow(merged) > 0)
  expect_true(all(rownames(merged) %in% cell_annotation_file$name))
})

test_that("data_prepare works with a precomputed nnPCA multi-score data set", {
  data(cell_annotation_file)
  data(nnPCA_Result_multiple)

  merged <- data_prepare(cell_annotation_file, nnPCA_Result_multiple,
    merge_colname = "name"
  )
  expect_s3_class(merged, "data.frame")
  # score columns from nnPCA_Result_multiple survive the merge
  expect_true(all(colnames(nnPCA_Result_multiple) %in% colnames(merged)))
  expect_true(nrow(merged) > 0)
})

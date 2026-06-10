test_that("correlate_sample_scores returns expected structure", {
  set.seed(42)
  ids <- paste0("s", seq_len(30))
  m1 <- data.frame(P1 = rnorm(30), P2 = rnorm(30), row.names = ids)
  m2 <- data.frame(EMT = rnorm(30), row.names = ids)

  res <- correlate_sample_scores(m1, m2, method = "pearson")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Correlation", "P_value") %in% colnames(res)))
  expect_equal(nrow(res), 2L)
})

test_that("correlate_sample_scores errors when no common samples", {
  m1 <- data.frame(P1 = rnorm(5), row.names = paste0("a", seq_len(5)))
  m2 <- data.frame(E = rnorm(5), row.names = paste0("b", seq_len(5)))
  expect_error(correlate_sample_scores(m1, m2))
})

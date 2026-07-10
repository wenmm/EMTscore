# Tests for the GMT utility helpers not already covered by test-gmt.R.
# (test-gmt.R covers the basic write_gmt/read_gmt round-trip and the unnamed
#  list error; here we extend to descriptions, multi-set files, and
#  filter_gmt_by_reference.)

test_that("read_gmt on the bundled test.gmt returns a gene data.frame", {
  skip_if_not_installed("GSEABase")
  gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

  genes <- read_gmt(gmt_file)
  expect_s3_class(genes, "data.frame")
  expect_equal(colnames(genes), "gene")
  expect_type(genes$gene, "character")
  expect_gt(nrow(genes), 0L)
  # genes are unique
  expect_equal(anyDuplicated(genes$gene), 0L)
})

test_that("write_gmt/read_gmt round-trips a multi-set file and dedups genes", {
  skip_if_not_installed("GSEABase")
  gene_sets <- list(
    SetA = c("TP53", "BRCA1", "EGFR"),
    # SetB shares EGFR with SetA -> read_gmt should return the union (unique)
    SetB = c("EGFR", "MYC", "VIM")
  )
  f <- tempfile(fileext = ".gmt")
  write_gmt(gene_sets, f)
  expect_true(file.exists(f))

  genes <- read_gmt(f)
  expect_setequal(genes$gene, c("TP53", "BRCA1", "EGFR", "MYC", "VIM"))
  expect_equal(nrow(genes), 5L)
})

test_that("write_gmt rejects a gene set that is not a character vector", {
  expect_error(
    write_gmt(list(BadSet = c(1, 2, 3)), tempfile(fileext = ".gmt")),
    "character vector"
  )
})

test_that("filter_gmt_by_reference writes a file and returns a stats data.frame", {
  skip_if_not_installed("GSEABase")
  ref_gmt <- system.file("extdata", "TianLab_collected_EMT_signatures.gmt",
    package = "EMTscore"
  )
  target_gmt <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt",
    package = "EMTscore"
  )
  output_gmt <- tempfile(fileext = ".gmt")

  stats <- suppressWarnings(suppressMessages(
    filter_gmt_by_reference(ref_gmt, target_gmt, output_gmt,
      cutoff = 0.5, keep_low_overlap = TRUE
    )
  ))

  expect_true(file.exists(output_gmt))
  expect_s3_class(stats, "data.frame")
  expect_equal(colnames(stats), c("name", "size", "overlap", "fraction"))
  expect_gt(nrow(stats), 0L)
  # fraction = overlap / size, and lies in [0, 1]
  expect_true(all(stats$fraction >= 0 & stats$fraction <= 1))
  expect_equal(stats$fraction, stats$overlap / stats$size)
})

test_that("filter_gmt_by_reference validates its numeric/logical arguments", {
  skip_if_not_installed("GSEABase")
  ref_gmt <- system.file("extdata", "TianLab_collected_EMT_signatures.gmt",
    package = "EMTscore"
  )
  target_gmt <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt",
    package = "EMTscore"
  )
  out <- tempfile(fileext = ".gmt")

  # cutoff outside [0, 1]
  expect_error(filter_gmt_by_reference(ref_gmt, target_gmt, out, cutoff = 2))
  # non-logical keep_low_overlap
  expect_error(
    filter_gmt_by_reference(ref_gmt, target_gmt, out, keep_low_overlap = "yes")
  )
})

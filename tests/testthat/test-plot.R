# Tests for the non-Seurat plotting helpers. These depend on Bioc/CRAN
# graphics packages, so each is guarded with skip_if_not_installed().

test_that("plot_heatmap_function returns a ComplexHeatmap Heatmap object", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  skip_if_not_installed("nsprcomp")
  data(geneExp)
  data(Tan_et_al_cell_line_M_signature)

  p <- plot_heatmap_function(t(geneExp), Tan_et_al_cell_line_M_signature)
  expect_s4_class(p, "Heatmap")
})

test_that("Execute_E_M_plot returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  data(data_for_plot)

  p <- Execute_E_M_plot(
    data_for_plot,
    E_colname = "Panchy_et_al_E_signature",
    M_colname = "Panchy_et_al_M_signature",
    celltype_colname = "celltype_annotation"
  )
  expect_s3_class(p, "ggplot")
})

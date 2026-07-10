#' EMTscore: Calculate 'EMT' Scores Based on 'Omics' Data
#'
#' The EMTscore package provides a toolbox to compute epithelial-mesenchymal
#' transition (EMT) scores from bulk and single-cell omics data using a range
#' of gene-set scoring algorithms and curated EMT gene sets, together with
#' visualization helpers for publication-quality figures.
#'
#' @keywords internal
#'
#' @importFrom Seurat AddModuleScore AddMetaData GetAssayData UpdateSeuratObject as.Seurat
#' @importFrom dplyr filter group_by summarise n
#' @import mclust
#' @import Matrix
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom stats cor sd complete.cases cor.test kmeans
#' @importFrom utils read.delim download.file
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom circlize colorRamp2
"_PACKAGE"

# Non-standard-evaluation variables used inside dplyr/foreach pipelines.
utils::globalVariables(c(
  "k", "mE", "mM", "sdE", "sdM", "mM1", "mM2", "sdM1", "sdM2"
))

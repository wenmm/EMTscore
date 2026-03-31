#' Compute EMT scores for one or multiple bulk expression matrices
#'
#' Calculate EMT (epithelial–mesenchymal transition) scores for bulk RNA-seq
#' expression data. The function accepts a single expression matrix (rows = genes,
#' columns = samples) or a named list of such matrices, and computes EMT scores
#' for each sample using one of the supported scoring methods.
#'
#' Supported methods:
#' \itemize{
#'   \item \code{"nnPCA"} — non-negative sparse PCA (via \code{nsprcomp})
#'   \item \code{"AUCell"} — AUCell AUC scoring
#'   \item \code{"ssGSEA"} — single-sample GSEA (GSVA)
#'   \item \code{"SCSE"} — Single-Cell Signature Explorer (adapted)
#'   \item \code{"JASMINE"} — JASMINE scoring
#' }
#'
#' @param expr_mat_list A matrix/data.frame (rows = genes, cols = samples) **or**
#'   a named list of such matrices. If a single matrix is provided it will be
#'   processed and a data.frame of scores returned. If a list is provided, each
#'   element will be processed and a named list of data.frames will be returned.
#' @param gmt_file Character. Path to a GMT file (gene sets). The GMT should
#'   contain EMT-related genes used by the scoring functions.
#' @param emt_name Character. Column name to use for the EMT score in the
#'   returned data.frames. Default: \code{"EMT_Score"}.
#' @param method Character. One of \code{"nnPCA"}, \code{"AUCell"},
#'   \code{"ssGSEA"}, \code{"SCSE"}, or \code{"JASMINE"}. Default selects the
#'   first value.
#' @param dimension Integer. Number of components for nnPCA (only used when
#'   \code{method = "nnPCA"}). Default: 1.
#'
#' @return If a single matrix is provided, a \code{data.frame} with rownames =
#'   sample names and one column named after \code{emt_name} containing EMT
#'   scores. If a list of matrices is provided, a named \code{list} of such
#'   \code{data.frame}s (one per input matrix) is returned.
#'
#' @details
#' - The function expects gene names as \code{rownames(mat)} and sample names as
#'   \code{colnames(mat)}. It will coerce inputs to numeric matrices internally.
#' - If none of the GMT genes are present in a matrix a warning is issued and
#'   that matrix returns \code{NULL} in the result list.
#' - The scoring helpers (\code{Execute_nnPCA}, \code{Execute_AUCell},
#'   \code{Execute_ssGSVA}, \code{Execute_SCSE}, \code{Execute_JAS}) are assumed
#'   available in the package namespace (or the environment) and must return a
#'   data.frame with sample rows and the score column name given in
#'   \code{score_names}.
#'
#' @examples
#' # single matrix
#' url <- "https://zenodo.org/records/18168504/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' res <- add_EMT_score_Bulk(expr_mat = geneExp, gmt_file, method = "AUCell")
#'
#' # multiple matrices (named list)
#' exprs <- list(TCGA = geneExp, GTEX = geneExp)
#' res_list <- add_EMT_score_Bulk(expr_mat_list = exprs, gmt_file, method = "nnPCA", dimension = 1)
#'
#' @export

add_EMT_score_Bulk <- function(expr_mat_list, gmt_file = NULL, 
                               emt_name, 
                               method = c("nnPCA", "AUCell", "ssGSEA", "SCSE", "JASMINE"), 
                               dimension) {
  method <- match.arg(method)
  
  if (is.null(gmt_file)) stop("gmt_file must be provided.")
  Genesets <- read_gmt(gmt_file)
  emt_genes <- Genesets$gene
  
  # 用来存放每个矩阵的结果
  results_list <- list()
  
  for (name in names(expr_mat_list)) {
    mat <- expr_mat_list[[name]]
    result_df <- data.frame(row.names = colnames(mat))
    
    if (method == "nnPCA") {
      scores <- Execute_nnPCA(mat, gmt_file, dimension, score_names = emt_name)
      result_df[[emt_name]] <- scores[rownames(result_df), emt_name]
      
    } else if (method == "AUCell") {
      scores <- Execute_AUCell(mat, gmt_file, score_names = emt_name)
      result_df[[emt_name]] <- scores[rownames(result_df), emt_name]
      
    } else if (method == "ssGSEA") {
      scores <- Execute_ssGSVA(mat, gmt_file, score_names = emt_name)
      result_df[[emt_name]] <- scores[rownames(result_df), emt_name]
      
    } else if (method == "SCSE") {
      scores <- Execute_SCSE(mat, gmt_file, score_names = emt_name)
      result_df[[emt_name]] <- scores[rownames(result_df), emt_name]
      
    } else if (method == "JASMINE") {
      scores <- Execute_JAS(mat, gmt_file, score_names = emt_name)
      result_df[[emt_name]] <- scores[rownames(result_df), emt_name]
    }
    
    results_list[[name]] <- result_df
  }
  
  return(results_list)
}



#' Execute Non-Negative Sparse PCA (nnPCA) Scoring
#'
#' Execute nnPCA Scoring from GMT File
#'
#' Read a GMT file and compute EMT scores for a gene expression matrix.
#'
#' @param exprMatrix Numeric matrix of gene expression, row names are genes, column names are samples.
#' @param Genesets Character string. Path to GMT file.
#' @param dimension An integer specifying the number of principal components to compute (default is 1).
#' @param score_names Character string. Column name for the scores.
#' @return Data frame of samples with nnPCA scores.
#'
#' @export
#' @import GSA stringr nsprcomp
#' 
#' @examples
#' url <- "https://zenodo.org/records/18168504/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' result <- Execute_nnPCA(geneExp, gmt_file,dimension = 1, score_names = "score")

Execute_nnPCA <- function(exprMatrix, Genesets, dimension, score_names)
{
  exprMatrix <- t(exprMatrix)
  genes <- unlist(read_gmt(Genesets)$gene)
  geneExp <- exprMatrix[, colnames(exprMatrix) %in% genes]
  
  if(ncol(geneExp) == 0){
    stop("No matching genes found in the expression matrix. Please check that the species and the capitalization of gene names in your expression matrix match those in the provided signature gene list.")
  }
  
  #pc_feature <- nsprcomp(as.matrix(geneExp), nneg=TRUE, ncomp=dimension, scale. = TRUE)
  pc_feature <- nsprcomp(as.matrix(geneExp), nneg=TRUE, ncomp=dimension)
  dfpc <- data.frame(pc_feature$x)
  dfpc[is.na(dfpc)] <- 0
  colnames(dfpc) <- score_names
  return(dfpc)
}

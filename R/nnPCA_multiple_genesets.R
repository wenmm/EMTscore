#' Compute nnPCA Scores for Multiple Gene Sets in Parallel
#'
#' This function reads a GMT gene set file and performs non-negative sparse PCA (nnPCA)
#' for each gene set in parallel, returning a combined data frame of pathway scores.
#'
#' @param exprMatrix Numeric matrix of gene expression, where rows are genes and columns are samples.
#' @param Genesets Character string. Path to the GMT file containing gene sets.
#' @param dimension Integer. Number of nnPCA components to compute.
#' @param cores Integer. Number of cores to use for parallel processing.
#'
#' @return A data frame containing the first nnPCA component scores for each pathway.
#' @export
#' @import foreach doParallel nsprcomp GSA
#'
#' @examples
#' url <- "https://zenodo.org/records/19487376/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
#' result <- Execute_nnPCA_parallel(geneExp, gmt_file, dimension = 1, cores = 1)

Execute_nnPCA_parallel <- function(exprMatrix, Genesets, dimension, cores) {
  # Load gene sets
  Genesets_obj <- GSA.read.gmt(Genesets)
  GSsize <- length(Genesets_obj$genesets)
  
  # Register parallel backend
  registerDoParallel(cores)
  
  # Parallel nnPCA computation
  CombinennPCA <- foreach(k = seq_len(GSsize), .combine = rbind, .errorhandling = "remove") %dopar% {
    genes <- unlist(Genesets_obj$genesets[k])
    pathwayName <- Genesets_obj$geneset.names[k]
    
    sub_expr <- exprMatrix[rownames(exprMatrix) %in% genes, , drop = FALSE]
    if (nrow(sub_expr) < 2) return(NULL)  # Skip small sets
    
    nnPCA_model <- nsprcomp(t(sub_expr), nneg = TRUE, ncomp = dimension)
    scores <- nnPCA_model$x[, 1]
    df <- data.frame(t(scores), stringsAsFactors = FALSE)
    df$Pathway <- pathwayName
    rownames(df) <- df$Pathway
    df$Pathway <- NULL
    df
  }
  
  return(t(CombinennPCA))
}

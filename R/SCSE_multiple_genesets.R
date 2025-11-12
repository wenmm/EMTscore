#' Compute SCSE Scores for Multiple Gene Sets
#'
#' This function computes SCSE (Single Cell Signature Explorer) scores for each gene set
#' in a GMT file using a given gene expression matrix (genes as rows, samples as columns).
#' It runs in parallel across multiple cores and returns a combined data frame
#' of SCSE scores for all pathways.
#'
#' @param exprMatrix Numeric matrix of gene expression, where rows are genes and columns are samples.
#' @param gmt_file Character string. Path to the GMT file containing gene sets.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing SCSE scores for each pathway across samples.
#' @export
#' @import foreach doParallel GSA stringr
#'
#' @examples
#' url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
#' result <- Execute_SCSE_parallel(data, gmt_file, cores = 2)

Execute_SCSE_parallel <- function(exprMatrix, gmt_file, cores) {
  stringsAsFactors <- FALSE

  registerDoParallel(cores)
  
  #---------------- SCSE function ----------------#
  SingleCellSigExplorer <- function(data, genes) {
    DataRanks <- data[rownames(data) %in% genes, , drop = FALSE]
    if (nrow(DataRanks) == 0) return(NULL)
    
    CumSum <- data.frame(colSums(DataRanks))
    colnames(CumSum)[1] <- "RawRankSum"
    CumSum$SampleID <- rownames(CumSum)
    
    TotalUMICount <- data.frame(colSums(data))
    colnames(TotalUMICount)[1] <- "TotalUMISum"
    TotalUMICount$SampleID <- rownames(TotalUMICount)
    
    FinalScore <- merge(CumSum, TotalUMICount, by = "SampleID")
    FinalScore$scSigExp <- FinalScore$RawRankSum / FinalScore$TotalUMISum * 100
    FinalScore[, c("SampleID", "scSigExp")]
  }
  
  #---------------- Wrapper to execute SCSE ----------------#
  Execute_SCSE <- function(data, Genesets, k) {
    genes <- unlist(Genesets$genesets[k])
    pathwayName <- Genesets$geneset.names[k]
    SCSE <- SingleCellSigExplorer(data, genes)
    if (is.null(SCSE)) return(NULL)
    
    SCSE_t <- as.data.frame(t(SCSE$scSigExp))
    names(SCSE_t) <- SCSE$SampleID
    SCSE_t$Pathway <- pathwayName
    SCSE_t
  }
  
  #---------------- Read gene sets ----------------#
  Genesets <- GSA.read.gmt(gmt_file)
  GSsize <- length(Genesets$genesets)
  
  #---------------- Parallel SCSE computation ----------------#
  Combine_SCSE <- foreach(k = 1:GSsize, .combine = rbind, .errorhandling = "remove") %dopar% {
    message("Processing gene set ", k, "/", GSsize, ": ", Genesets$geneset.names[k])
    Execute_SCSE(exprMatrix, Genesets, k)
  }
  rownames(Combine_SCSE) <- Combine_SCSE$Pathway
  Combine_SCSE$Pathway <- NULL
  return(t(Combine_SCSE))
}

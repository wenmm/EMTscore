#' Compute AUCell Scores for Multiple Gene Sets
#'
#' This function computes AUCell scores for each gene set in a GMT file
#' using a given gene expression matrix (genes as rows, samples as columns).
#' It runs in parallel across multiple cores and returns a combined data frame
#' of AUCell scores for all pathways.
#'
#' @param exprMatrix Numeric matrix or data frame of gene expression, where row names are genes and column names are samples.
#' @param gmt_file Character string. Path to the GMT file containing gene sets.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing AUCell scores for each pathway across samples.
#' @export
#' @import AUCell foreach doParallel GSA stringr
#'
#' @examples
#' url <- "https://zenodo.org/records/19487376/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
#' result <- Execute_AUCell_parallel(geneExp, gmt_file, cores = 2)

Execute_AUCell_parallel <- function(exprMatrix, gmt_file, cores) {
  stringsAsFactors <- FALSE
  registerDoParallel(cores)
  
  #---------------- AUCell scoring function for one gene set ----------------#
  AUCellfunc <- function(exprMatrix, genes) {
    out <- tryCatch({
      cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
      geneSets <- list(geneSet1 = genes)
      cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
      cellsAUCretrieve <- getAUC(cells_AUC)
      cells_AUCellScores <- t(cellsAUCretrieve)
      cells_AUCellScores <- data.frame(rownames(cells_AUCellScores), cells_AUCellScores)
      colnames(cells_AUCellScores) <- c("SampleID", "AUCell")
      row.names(cells_AUCellScores) <- NULL
      return(cells_AUCellScores)
    },
    error = function(cond) {
      cells_AUCellScores <- data.frame(SampleID = "NA", AUCell = "NA")
      return(cells_AUCellScores)
    })
    return(out)
  }
  
  #---------------- Function to execute AUCell ----------------#
  Execute_AUCell <- function(exprMatrix, Genesets, k) {
    genes <- unlist(Genesets$genesets[k])
    pathwayName <- Genesets$geneset.names[k]
    AUCellMethod <- AUCellfunc(exprMatrix, genes)
    if (nrow(AUCellMethod) == ncol(exprMatrix)) {
      Samples <- AUCellMethod$SampleID
      AUCellMethod <- data.frame(AUCellMethod[, -1])
      AUCellMethod <- data.frame(t(AUCellMethod))
      names(AUCellMethod) <- Samples
      AUCellMethod$Pathway <- pathwayName
      return(AUCellMethod)
    }
  }
  
  #---------------- Read gene sets ----------------#
  Genesets <- GSA.read.gmt(gmt_file)
  GSsize <- length(Genesets$genesets)
  
  #---------------- Run AUCell in parallel ----------------#
  Combine_AUCellResult <- foreach(k = seq_len(GSsize), .combine = rbind, .errorhandling = "remove") %dopar% {
    message("Processing pathway ", k, "/", GSsize, ": ", Genesets$geneset.names[k])
    Execute_AUCell(exprMatrix, Genesets, k)
  }
  
  rownames(Combine_AUCellResult) <- Combine_AUCellResult$Pathway
  Combine_AUCellResult$Pathway <- NULL
  Combine_AUCellResult
  return(t(Combine_AUCellResult))
}

#' Compute JASMINE Scores for Multiple Gene Sets in Parallel
#'
#' This function calculates JASMINE scores for gene sets defined in a GMT file.
#' Each pathway score combines normalized rank mean and odds ratio of gene expression across samples.
#' Parallel processing is used to speed up computation across gene sets.
#'
#' @param exprMatrix Numeric matrix of gene expression, rows are genes, columns are samples.
#' @param gmt_file Character string. Path to the GMT file containing gene sets.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing JASMINE scores for each pathway across samples.
#' @export
#' @importFrom GSA GSA.read.gmt
#'
#' @examples
#' data(geneExp)
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' result <- Execute_JASMINE_parallel(geneExp, gmt_file, cores = 1)
Execute_JASMINE_parallel <- function(exprMatrix, gmt_file, cores) {

  #----------- Internal helper functions -----------#
  RankCalculation <- function(x, genes) {
    subdata <- x[x != 0]
    if (length(subdata) == 0) {
      return(0)
    }
    DataRanksUpdated <- rank(subdata)
    DataRanksSigGenes <- DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]
    CumSum <- if (length(DataRanksSigGenes) == 0) 0 else mean(DataRanksSigGenes, na.rm = TRUE)
    FinalRawRank <- CumSum / length(subdata)
    return(FinalRawRank)
  }

  ORCalculation <- function(data, genes) {
    GE <- data[rownames(data) %in% genes, , drop = FALSE]
    NGE <- data[!rownames(data) %in% genes, , drop = FALSE]
    SigGenesExp <- apply(GE, 2, function(x) sum(x != 0))
    NSigGenesExp <- apply(NGE, 2, function(x) sum(x != 0))
    SigGenesNE <- nrow(GE) - SigGenesExp
    SigGenesNE[SigGenesNE == 0] <- 1
    NSigGenesExp[NSigGenesExp == 0] <- 1
    NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp) - SigGenesNE
    denom <- SigGenesNE * NSigGenesExp
    denom[denom == 0] <- 1
    OR <- (SigGenesExp * NSigGenesNE) / denom
    OR[is.na(OR) | is.infinite(OR)] <- 0
    return(OR)
  }

  NormalizationJAS <- function(x) {
    if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) {
      return(rep(0, length(x)))
    }
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }

  JASMINE_single <- function(data, genes, pathwayName) {
    idx <- match(genes, rownames(data))
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) {
      return(NULL)
    }
    RM <- apply(data, 2, function(x) RankCalculation(x, genes))
    OR <- ORCalculation(data, genes)
    RM <- NormalizationJAS(RM)
    OR <- NormalizationJAS(OR)
    JAS_Scores <- (RM + OR) / 2
    FinalScores <- data.frame(t(JAS_Scores))
    colnames(FinalScores) <- colnames(data)
    FinalScores$Pathway <- pathwayName
    return(FinalScores)
  }

  #----------- Read gene sets -----------#
  Genesets <- GSA.read.gmt(gmt_file)
  GSsize <- length(Genesets$genesets)

  #----------- Parallel computation -----------#
  Combine_JASMINE <- bind_genesets_parallel(GSsize, function(k) {
    genes <- unlist(Genesets$genesets[k])
    pathwayName <- Genesets$geneset.names[k]
    res <- JASMINE_single(exprMatrix, genes, pathwayName)
    res
  }, cores)


  # Remove empty rows
  if (nrow(Combine_JASMINE) == 0) {
    message("NULL result")
  } else {
    Combine_JASMINE <- Combine_JASMINE[rowSums(is.na(Combine_JASMINE)) != ncol(Combine_JASMINE), ]
  }

  rownames(Combine_JASMINE) <- Combine_JASMINE$Pathway
  Combine_JASMINE$Pathway <- NULL
  return(t(Combine_JASMINE))
}

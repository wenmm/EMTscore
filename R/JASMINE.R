### Rank Calculation function
RankCalculation <- function(x, genes) {
  # Keep only expressed genes
  subdata <- x[x != 0]
  
  if (length(subdata) == 0) {
    return(0)  # no expression in this cell
  }
  
  DataRanksUpdated <- rank(subdata)
  DataRanksSigGenes <- DataRanksUpdated[which(names(DataRanksUpdated) %in% genes)]
  
  if (length(DataRanksSigGenes) == 0) {
    CumSum <- 0
  } else {
    CumSum <- mean(DataRanksSigGenes, na.rm = TRUE)
  }
  
  FinalRawRank <- CumSum / length(subdata)
  return(FinalRawRank)
}



ORCalculation <- function(data, genes) {
  GE <- data[which(rownames(data) %in% genes), , drop = FALSE]
  NGE <- data[-which(rownames(data) %in% genes), , drop = FALSE]
  
  SigGenesExp <- apply(GE, 2, function(x) length(x[x != 0]))
  NSigGenesExp <- apply(NGE, 2, function(x) length(x[x != 0]))
  
  SigGenesNE <- nrow(GE) - SigGenesExp
  SigGenesNE[SigGenesNE == 0] <- 1
  
  NSigGenesExp[NSigGenesExp == 0] <- 1
  NSigGenesNE <- nrow(data) - (NSigGenesExp + SigGenesExp)
  NSigGenesNE <- NSigGenesNE - SigGenesNE
  
  denom <- SigGenesNE * NSigGenesExp
  denom[denom == 0] <- 1  # avoid divide-by-zero
  
  OR <- (SigGenesExp * NSigGenesNE) / denom
  OR[is.na(OR) | is.infinite(OR)] <- 0
  return(OR)
}


NormalizationJAS <- function(JAS_Scores) {
  if (max(JAS_Scores, na.rm = TRUE) == min(JAS_Scores, na.rm = TRUE)) {
    return(rep(0, length(JAS_Scores)))
  } else {
    return((JAS_Scores - min(JAS_Scores, na.rm = TRUE)) /
             (max(JAS_Scores, na.rm = TRUE) - min(JAS_Scores, na.rm = TRUE)))
  }
}


JASMINE <- function(data, genes) {
  idx = match(genes, rownames(data))
  idx = idx[!is.na(idx)]
  if (length(idx) > 1) {
    RM = apply(data, 2, function(x) RankCalculation(x, genes))
    OR = ORCalculation(data, genes)
    RM = NormalizationJAS(RM)
    OR = NormalizationJAS(OR)
    JAS_Scores = (RM + OR) / 2
    FinalScores = data.frame(names(RM), JAS_Scores)
    colnames(FinalScores)[1] = "SampleID"
    Samples = FinalScores$SampleID
    FinalScores = FinalScores[, -1]
    FinalScores = data.frame(t(FinalScores))
    names(FinalScores) = as.character(Samples)
    
    return(FinalScores)
  }
}

#' Execute JASMINE Scoring from GMT File
#'
#' Read a GMT file and compute JASMINE scores for a gene expression matrix.
#'
#' @param exprMatrix Numeric matrix of gene expression, row names are genes, column names are samples.
#' @param Genesets Character string. Path to GMT file.
#' @param score_names Character string. Column name for the scores.
#' @return Data frame of samples with JASMINE scores.
#'
#' @export
#' @import GSA stringr
#' 
#' @examples
#' url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' result <- Execute_JAS(geneExp, gmt_file, score_names = "score")
Execute_JAS <- function(exprMatrix, Genesets, score_names){
  genes <- unlist(read_gmt(Genesets)$gene)
  subsetJAS <- JASMINE(exprMatrix, genes)
  subsetJAS = data.frame(t(subsetJAS))
  names(subsetJAS)= score_names
  return(subsetJAS)
}
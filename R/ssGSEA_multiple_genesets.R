#' Compute ssGSEA Scores for Multiple Gene Sets in Parallel
#'
#' This function calculates ssGSEA scores for gene sets defined in a GMT file
#' using the ssGSEA package. Parallel processing is used to speed up computation across pathways.
#'
#' @param exprMatrix Numeric matrix of gene expression, rows are genes, columns are samples.
#' @param gmt_file Character string. Path to the GMT file containing gene sets.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing ssGSEA scores for each pathway across samples. Each row corresponds to a pathway with an additional `Pathway` column.
#' @export
#' @importFrom GSVA ssgseaParam gsva
#' @importFrom GSA GSA.read.gmt
#'
#' @examples
#' data(geneExp)
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' result <- Execute_ssGSEA_parallel(geneExp, gmt_file, cores = 1)
Execute_ssGSEA_parallel <- function(exprMatrix, gmt_file, cores) {

  #----------- Read gene sets -----------#
  Genesets <- GSA.read.gmt(gmt_file)
  GSsize <- length(Genesets$genesets)

  #----------- Wrapper for one gene set -----------#
  Compute_ssGSEA <- function(exprMatrix, Genesets, k) {
    geneslist <- unlist(Genesets$genesets[k])
    pathwayName <- Genesets$geneset.names[k]

    gsva_geneSets <- GSEABase::GeneSetCollection(
      GSEABase::GeneSet(geneslist, setName = pathwayName)
    )
    gsvaPar <- GSVA::ssgseaParam(exprMatrix, geneSets = gsva_geneSets)
    Result <- gsva(gsvaPar, verbose = FALSE)
    Result <- data.frame(Result)
    Result$Pathway <- pathwayName
    return(Result)
  }

  #----------- Parallel computation -----------#
  Combine_ssGSEA <- bind_genesets_parallel(GSsize, function(k) {
    Compute_ssGSEA(exprMatrix, Genesets, k)
  }, cores)

  # Remove empty rows if any
  if (nrow(Combine_ssGSEA) == 0) {
    message("NULL result")
  } else {
    Combine_ssGSEA <- Combine_ssGSEA[rowSums(is.na(Combine_ssGSEA)) != ncol(Combine_ssGSEA), ]
  }
  rownames(Combine_ssGSEA) <- Combine_ssGSEA$Pathway
  Combine_ssGSEA$Pathway <- NULL
  return(t(Combine_ssGSEA))
}

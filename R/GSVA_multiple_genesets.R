#' Compute GSVA Scores for Multiple Gene Sets in Parallel
#'
#' This function calculates GSVA scores for gene sets defined in a GMT file
#' using the GSVA package. Parallel processing is used to speed up computation across pathways.
#'
#' @param exprMatrix Numeric matrix of gene expression, rows are genes, columns are samples.
#' @param gmt_file Character string. Path to the GMT file containing gene sets.
#' @param cores Integer. Number of CPU cores to use for parallel processing.
#'
#' @return A data frame containing GSVA scores for each pathway across samples. Each row corresponds to a pathway with an additional `Pathway` column.
#' @export
#' @import GSVA GSA foreach doParallel stringr
#'
#' @examples
#' url <- "https://zenodo.org/records/19487376/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
#' result <- Execute_GSVA_parallel(geneExp, gmt_file, cores = 2)

Execute_GSVA_parallel <- function(exprMatrix, gmt_file, cores) {
  
  registerDoParallel(cores)
  
  #----------- Read gene sets -----------#
  Genesets <- GSA.read.gmt(gmt_file)
  GSsize <- length(Genesets$genesets)
  
  #----------- Wrapper for one gene set -----------#
  Compute_GSVA <- function(exprMatrix, Genesets, k) {
    geneslist <- unlist(Genesets$genesets[k])
    pathwayName <- Genesets$geneset.names[k]
    
    gsva_geneSets <- GSEABase::GeneSetCollection(
      GSEABase::GeneSet(geneslist, setName = pathwayName)
    )
    gsvaPar <- gsvaParam(exprMatrix, geneSets = gsva_geneSets)
    res <- tryCatch({
      Result <- gsva(gsvaPar, verbose = FALSE)
      Result <- data.frame(Result)
      Result$Pathway <- pathwayName
      return(Result)
    }, error = function(e) {
      NULL
    })
    return(res)
  }
  
  #----------- Parallel computation -----------#
  Combine_GSVA <- foreach(k = seq_len(GSsize), .combine = rbind, .errorhandling = "remove") %dopar% {
    Compute_GSVA(exprMatrix, Genesets, k)
  }
  
  # Remove empty rows if any
  if (nrow(Combine_GSVA) == 0) {
    message("NULL result")
  } else {
    Combine_GSVA <- Combine_GSVA[rowSums(is.na(Combine_GSVA)) != ncol(Combine_GSVA), ]
  }
  rownames(Combine_GSVA) <- Combine_GSVA$Pathway
  Combine_GSVA$Pathway <- NULL
  return(t(Combine_GSVA))
}

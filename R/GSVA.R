#' Execute GSVA Scoring from GMT File
#'
#' Read a GMT file and compute EMT scores for a gene expression matrix.
#'
#' @param exprMatrix Numeric matrix of gene expression, row names are genes, column names are samples.
#' @param Genesets Character string. Path to GMT file.
#' @param score_names Character string. Column name for the scores.
#' @return Data frame of samples with GSVA scores.
#'
#' @export
#' @import GSA stringr GSVA
#' 
#' @examples
#' url <- "https://zenodo.org/records/18168504/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' result <- Execute_GSVA(geneExp, gmt_file, score_names = "score")

Execute_GSVA <- function(exprMatrix,Genesets,score_names){
geneslist = unlist(read_gmt(Genesets)$gene)
gsva_geneSets <- GSEABase::GeneSetCollection(
  GSEABase::GeneSet(geneslist, setName = score_names)
)
gsvaPar <- gsvaParam(exprMatrix, geneSets = gsva_geneSets)
Result = gsva(gsvaPar, verbose=FALSE)     
Result = data.frame(t(Result))
colnames(Result) <- score_names
return(Result)
}

########################### AUCell Function to compute the scores ###############################################

## Function Input: single-cell gene expression data and genes in one gene sets

## Function Output: Scores for each cell/sample ID


AUCellfunc <- function(exprMatrix,genes){

out <- tryCatch(
{
cells_rankings <- AUCell_buildRankings(exprMatrix,plotStats=FALSE)
geneSets <- list(geneSet1=genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cellsAUCretrieve <- getAUC(cells_AUC)
cells_AUCellScores <- t(cellsAUCretrieve)
cells_AUCellScores <- data.frame(rownames(cells_AUCellScores),cells_AUCellScores)
colnames(cells_AUCellScores) <- c('SampleID','AUCell')
row.names(cells_AUCellScores) <- NULL
return(cells_AUCellScores)

},

error = function(cond) {
cells_AUCellScores <- data.frame(SampleID = "NA", AUCell = "NA")
return(cells_AUCellScores)
}
)
return(out)
}


#' Execute AUCell Scoring from GMT File
#'
#' Read a GMT file and compute EMT scores for a gene expression matrix.
#'
#' @param exprMatrix Numeric matrix of gene expression, row names are genes, column names are samples.
#' @param Genesets Character string. Path to GMT file.
#' @param score_names Character string. Column name for the scores.
#' @return Data frame of samples with AUCell scores.
#'
#' @export
#' @import GSA stringr AUCell
#' 
#' @examples
#' url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' result <- Execute_AUCell(geneExp, gmt_file, score_names = "score")
Execute_AUCell <- function(exprMatrix,Genesets,score_names)
{
	genes = unlist(read_gmt(Genesets)$gene)
	AUCellMethod = AUCellfunc(exprMatrix,genes)
	if(nrow(AUCellMethod) == ncol(exprMatrix))
	{
		Samples = AUCellMethod$SampleID
		AUCellMethod = data.frame(AUCellMethod[,-1])
		AUCellMethod = data.frame(t(AUCellMethod))
		names(AUCellMethod) = Samples
		result = data.frame(t(AUCellMethod))
		colnames(result) <- score_names
		return(result)
	}
}

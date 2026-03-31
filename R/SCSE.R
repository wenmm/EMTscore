################ SCSE function ############ 

## Function Input: single-cell gene expression data and genes in one gene sets (This function is implemented using the equation described in Pont et al.,2019)

## Function Output: Scores for each cell/sample ID

SingleCellSigExplorer <- function(data,genes)
{
	DataRanks = data[which(rownames(data) %in% genes),]
	if(length(nrow(DataRanks))!=0)
	{
		CumSum = data.frame(colSums(DataRanks, na.rm = FALSE, dims = 1))
		colnames(CumSum)[1]='RawRankSum'
		SampleID = rownames(CumSum)
		CumSum1 = data.frame(SampleID,CumSum)
		row.names(CumSum1)=NULL
		TotalUMICount = data.frame(colSums(data, na.rm = FALSE, dims = 1))
		colnames(TotalUMICount)[1]='TotalUMISum'
		SampleID = rownames(TotalUMICount)
		TotalUMICount = data.frame(SampleID,TotalUMICount)
		row.names(TotalUMICount)=NULL
		FinalScore = merge(CumSum1, TotalUMICount, by='SampleID')
		FinalScore$SingleScorer = FinalScore$RawRankSum/FinalScore$TotalUMISum*100
		colnames(FinalScore)[4] = 'scSigExp'
		FinalScore = FinalScore[,c(1,4)]
		return(FinalScore)
	} 
}

#' Execute SCSE Scoring from GMT File
#'
#' Read a GMT file and compute EMT scores for a gene expression matrix.
#'
#' @param exprMatrix Numeric matrix of gene expression, row names are genes, column names are samples.
#' @param Genesets Character string. Path to GMT file.
#' @param score_names Character string. Column name for the scores.
#' @return Data frame of samples with SCSE scores.
#'
#' @export
#' @import GSA stringr
#' 
#' @examples
#' url <- "https://zenodo.org/records/18168504/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' result <- Execute_SCSE(geneExp, gmt_file, score_names = "score")

Execute_SCSE <- function(exprMatrix,Genesets,score_names)
{
	genes = unlist(read_gmt(Genesets)$gene)
	SCSE = SingleCellSigExplorer(exprMatrix,genes)
	if(length(SCSE)!=0)
	{
		Samples = SCSE$SampleID
		SCSE = data.frame(SCSE[,-1])
		SCSE = data.frame(t(SCSE))
		names(SCSE)= Samples
		SCSE = data.frame(t(SCSE))
		colnames(SCSE)= score_names
		return(SCSE)
	}
}

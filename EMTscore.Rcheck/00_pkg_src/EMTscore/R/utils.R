#' Read GMT File and Extract Genes
#'
#' This function reads a GMT (Gene Matrix Transposed) file and extracts all unique genes
#' across all gene sets.
#'
#' @param fname Character string. Path to the GMT file.
#'
#' @return A data frame with a single column \code{gene}, containing all unique genes.
#'
#' @examples
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#'
#' genes <- read_gmt(gmt_file)
#'
#' @export

read_gmt <- function(fname){
  gmt_lines <- readLines(fname)
  gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split="\t")))
  gmt_genes <- lapply(gmt_list, function(x) x[3:length(x)])
  genes <- unique(unlist(gmt_genes))
  return(data.frame(gene = genes, stringsAsFactors = FALSE))
}



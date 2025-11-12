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

#' Write a Gene Set Collection to a GMT File
#'
#' This function exports a list of gene sets into a GMT (Gene Matrix Transposed) file format.  
#' Each line of the GMT file corresponds to one gene set, with the structure:
#' By default, the gene set name is used as both the name and description fields.
#'
#' @param gene_sets A named list where each element is a character vector of gene identifiers.
#'   The list names correspond to gene set names.
#' @param file A character string specifying the path to the output GMT file.
#'
#' @details
#' The GMT format is widely used by tools such as GSEA and Enrichr for representing gene set collections.
#' Each line contains the set name, a description (often a URL or repeated name), and one or more gene identifiers separated by tabs.
#'
#' @return
#' This function writes a GMT file to disk and returns \code{NULL} invisibly.
#'
#' @examples
#' gene_sets <- list(
#'   E = c("CDH1", "VIM", "SNAI1", "ZEB1"),
#'   M = c("CDK1", "CCNB1", "CDC20")
#' )
#' write_gmt(gene_sets, file = "EM_signature.gmt")
#'
#' @export
write_gmt <- function(gene_sets, file) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
  }
  # Construct GMT lines: "name description gene1 gene2..."
  gmt_lines <- lapply(names(gene_sets), function(name) {
    genes <- gene_sets[[name]]
    if (!is.character(genes)) {
    }
    # Use name as both set name and description (second field)
    line <- paste(c(name, name, genes), collapse = "\t")
    return(line)
  })
  # Write to file
  writeLines(unlist(gmt_lines), con = file, useBytes = TRUE)
}


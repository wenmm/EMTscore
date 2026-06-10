#' Write Gene Sets to a GMT File
#'
#' Convert multiple gene lists into a single GMT (Gene Matrix Transposed) formatted file.  
#' Each element of the input list represents one gene set, where the element name  
#' is the gene set name and the element value is a vector of gene symbols.
#'
#' @param gene_sets A named list of gene sets. Each element should be a character vector  
#'        containing gene symbols, and each list name will be used as the gene set name.  
#'        Example: \code{list(Pathway1 = c("TP53", "BRCA1"), Pathway2 = c("EGFR", "MYC"))}
#' @param file A character string specifying the output GMT file path (e.g. "output.gmt").
#'
#' @return No return value. The function writes a GMT file to the specified path.
#' @export
#'
#' @examples
#' gene_sets <- list(
#'   EMT_High = c("VIM", "ZEB1", "SNAI1"),
#'   EMT_Low  = c("CDH1", "OCLN", "DSP")
#' )
#' write_gmt(gene_sets, file = tempfile(fileext = ".gmt"))
#'
write_gmt <- function(gene_sets, file) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a *named* list of gene vectors.")
  }
  
  # Construct GMT lines: "name<TAB>description<TAB>gene1<TAB>gene2..."
  gmt_lines <- lapply(names(gene_sets), function(name) {
    genes <- gene_sets[[name]]
    if (!is.character(genes)) {
      stop("Each gene set must be a character vector. Problem at: ", name)
    }
    # Use name as both set name and description (second field)
    line <- paste(c(name, name, genes), collapse = "\t")
    return(line)
  })
  
  # Write to file
  writeLines(unlist(gmt_lines), con = file, useBytes = TRUE)
}

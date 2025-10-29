#' Compute Signature Score from Bulk Matrix or Seurat Object
#'
#' This function computes signature scores using Spearman correlation 
#' between a predefined gene signature and expression data. It supports 
#' both bulk expression matrices and Seurat objects.
#'
#' @param data A bulk gene expression matrix (genes × samples) or a Seurat object.
#' @param signature_file Path to the signature weight file (two-column .tsv file).
#' @param row_name Character string. The name to assign to the resulting score column.
#'
#' @return A numeric matrix (for bulk input) or a Seurat object with a new metadata column (for Seurat input).
#' @examples
#' 
#' url <- https://zenodo.org/records/17438655/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
#' compute_Signature_score(geneExp, signature_file, "Stemness_Score")

compute_Signature_score <- function(data, signature_file, row_name) {
  
  if (!file.exists(signature_file)) stop("Signature file does not exist.")
  
  # Load signature weights
  w <- read.delim(signature_file, header = FALSE, row.names = 1) %>% as.matrix() %>% drop()
  
  if (inherits(data, "Seurat")) {
    # ----- Seurat object -----
    expr_mat <- Seurat::GetAssayData(data, assay = "RNA", slot = "data")
    w <- w[rownames(expr_mat)]
    
    s <- apply(expr_mat, 2, function(z) cor(z, w, method = "sp", use = "complete.obs"))
    s <- s - min(s)
    s <- s / max(s)
    result <- cbind(s)
    
    data@meta.data[[row_name]] <- result[rownames(data@meta.data)]
    message(paste("Added signature score to Seurat metadata as:", row_name))
    return(data)
    
  } else if (is.matrix(data) || is.data.frame(data)) {
    # ----- Bulk matrix -----
    expr_mat <- as.matrix(data)
    w <- w[rownames(expr_mat)]
    
    s <- apply(expr_mat, 2, function(z) cor(z, w, method = "sp", use = "complete.obs"))
    s <- s - min(s)
    s <- s / max(s)
    
    result <- cbind(s)
    colnames(result) <- row_name
    return(result)
    
  } else {
    stop("Input must be either a Seurat object or a numeric expression matrix.")
  }
}

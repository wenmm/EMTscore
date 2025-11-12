#' Compute Signature Score from Bulk Matrix or Seurat Object
#'
#' This function computes signature scores using Spearman correlation 
#' between a predefined gene signature and expression data. It supports 
#' both bulk expression matrices and Seurat objects.
#'
#' @param data A bulk gene expression matrix (row names are gene names, and the column names are sample) or a Seurat object.
#' @param signature_file Path to the signature weight file (two-column .tsv file).
#' @param score_name Character string. The name to assign to the resulting score column.
#'
#' @return A numeric matrix (for bulk input) or a Seurat object with a new metadata column (for Seurat input).
#' @examples 
#' url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
#' destfile <- tempfile(fileext = ".rda")
#' download.file(url, destfile, mode = "wb")
#' load(destfile)
#' signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
#' compute_Signature_score(geneExp, signature_file, "Stemness_Score")
#' @export

compute_Signature_score <- function(data, signature_file, score_name) {
  
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
    
    data@meta.data[[score_name]] <- result[rownames(data@meta.data)]
    return(data)
    
  } else if (is.matrix(data) || is.data.frame(data)) {
    # ----- Bulk matrix -----
    expr_mat <- as.matrix(data)
    w <- w[rownames(expr_mat)]
    
    s <- apply(expr_mat, 2, function(z) cor(z, w, method = "sp", use = "complete.obs"))
    s <- s - min(s)
    s <- s / max(s)
    
    result <- cbind(s)
    colnames(result) <- score_name
    return(result)
    
  } else {
    stop("Input must be either a Seurat object or a numeric expression matrix.")
  }
}


#' Compute Multiple Signature Scores for Single-Cell Seurat Objects
#'
#' This function calculates one or more gene signature scores for single-cell datasets 
#' stored as Seurat objects. Each score is computed using Spearman correlation 
#' between gene expression profiles and signature weights, followed by min-max normalization.
#' The resulting scores are added to the Seurat object metadata.
#'
#' @param seurat_files A character vector of paths to Seurat object `.rds` files.
#' Each file will be read, updated, and annotated with computed signature scores.
#' @param signature_files A character vector of paths to signature weight files 
#' (two-column `.tsv` or `.txt` files with gene symbols and weights). 
#' Multiple files can be provided.
#' @param score_names A character vector of names for the resulting score columns. 
#' Must have the same length as `signature_files`.
#'
#' @details 
#' Each signature score is calculated by computing the Spearman correlation 
#' between the normalized expression of genes and their corresponding signature weights. 
#' The resulting scores are then scaled to a 0-1 range and stored in the Seurat 
#' object metadata under the specified names.
#'
#' @return 
#' A named list of Seurat objects, each containing new metadata columns corresponding 
#' to the computed signature scores.
#'
#'
#' @importFrom Seurat GetAssayData UpdateSeuratObject
#' @importFrom stats cor
#' @importFrom utils read.delim
#' @importFrom magrittr %>%
#' @export
compute_Signature_score_SingleCell <- function(seurat_files, signature_files, score_names) {
  # Check input consistency
  if (length(signature_files) != length(score_names)) {
    stop("The number of signature files must match the number of score names.")
  }
  
  # Load all signature weight files
  weights_list <- lapply(signature_files, function(f) {
    w <- read.delim(f, header = FALSE, row.names = 1) %>% as.matrix() %>% drop()
    return(w)
  })
  
  # Helper function: compute normalized Spearman correlation score
  compute_score <- function(expr_mat, w) {
    w <- w[rownames(expr_mat)]
    s <- apply(expr_mat, 2, function(z) cor(z, w, method = "sp", use = "complete.obs"))
    s <- s - min(s)
    s <- s / max(s)
    return(s)
  }
  
  # Assign file names
  names(seurat_files) <- gsub("\\.rds$", "", basename(seurat_files))
  
  # Process each Seurat object
  seurat_list <- lapply(names(seurat_files), function(name) {
    obj <- readRDS(seurat_files[name])
    obj <- UpdateSeuratObject(obj)
    
    if (inherits(obj, "Seurat")) {
      expr_mat <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
      
      # Loop over all signatures
      for (i in seq_along(weights_list)) {
        score <- compute_score(expr_mat, weights_list[[i]])
        obj@meta.data[[score_names[i]]] <- score[rownames(obj@meta.data)]
      }
      return(obj)
    }
  })
  
  names(seurat_list) <- names(seurat_files)
  return(seurat_list)
}



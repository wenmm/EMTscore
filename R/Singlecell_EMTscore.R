#' Add EMT Scores to Seurat Objects
#'
#' This function calculates epithelial–mesenchymal transition (EMT) scores 
#' for one or more Seurat objects using a variety of scoring methods. 
#' Supported methods include \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"}, 
#' \code{"ssGSEA"}, \code{"SCSE"}, and \code{"JASMINE"}.
#'
#' @param seurat_files A character vector of file paths to Seurat objects (\code{.rds} files).
#' @param gmt_file Character string. Path to the gene set file in GMT format. 
#'   Must include an EMT gene signature.
#' @param emt_name Character string. Name of the EMT score column to be added 
#'   to the Seurat object metadata. Default: \code{"EMT_Score"}.
#' @param method Character string. Scoring method to use. 
#'   Must be one of: \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"}, 
#'   \code{"ssGSEA"}, \code{"SCSE"}, or \code{"JASMINE"}.
#' @param nnPCA_dim Integer. Dimension to use for nnPCA scoring (default = 1).
#'
#' @details
#' - \strong{Seurat}: Uses \code{Seurat::AddModuleScore} to compute module scores.  
#' - \strong{nnPCA}: Computes nonlinear PCA–based EMT scores.  
#' - \strong{AUCell}: Calculates the Area Under the Curve (AUC) for EMT gene set enrichment.  
#' - \strong{ssGSEA}: Runs single-sample GSEA for EMT signatures.  
#' - \strong{SCSE}: Uses the SCSE scoring method for EMT genes.  
#' - \strong{JASMINE}: Computes EMT scores via the JASMINE algorithm.  
#'
#' Each processed Seurat object will have a new column added to 
#' its \code{meta.data} slot containing EMT scores.
#'
#' @return A named list of Seurat objects with updated metadata, 
#'   each including an EMT score column.
#'
#' @examples
#' rds_file <- system.file("extdata", "example_seurat_obj.rds", package = "EMTscore")
#' seurat_files <- c("/home/hwen6/project/EMTscore/2025_09_24/EMTscore/inst/extdata/example_seurat_obj.rds", "/home/hwen6/project/EMTscore/2025_09_24/EMTscore/inst/extdata/example_seurat_obj.rds")
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#'
#' # Using nnPCA
#' seurat_list <- add_EMT_score(seurat_files, gmt_file, emt_name = "EMT_score", method = "nnPCA")
#'
#' @import Matrix patchwork dplyr slingshot velocyto.R ggplot2 Rmagic SeuratWrappers DoubletFinder monocle Seurat data.table nsprcomp
#' @export
#' 

add_EMT_score <- function(seurat_files, gmt_file = NULL,
                          emt_name, method = c("Seurat","nnPCA","AUCell","ssGSEA","SCSE","JASMINE"),
                          nnPCA_dim = 1) {
  
  method <- match.arg(method)
  
  # read geneset（用于 Seurat, nnPCA, AUCell, ssGSEA, SCSE）
  if (is.null(gmt_file)) stop("gmt_file must be provided.")
  Genesets <- read_gmt(gmt_file)
  
  # add label
  names(seurat_files) <- gsub(".rds", "", basename(seurat_files))
  
  seurat_list <- lapply(names(seurat_files), function(name){
    obj <- readRDS(seurat_files[name])
    obj <- UpdateSeuratObject(obj)
    emt_genes <- Genesets$gene
    
    if (method == "Seurat") {
      obj <- AddModuleScore(obj, features = list(emt_genes), name = "EMT_Score")
      colnames(obj@meta.data)[colnames(obj@meta.data) == "EMT_Score1"] <- emt_name
      
      if (!(emt_name %in% colnames(obj@meta.data))) stop(paste("Score column", emt_name, "not found."))
      
    } else if (method == "nnPCA") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      nnPCA_Result <- Execute_nnPCA(geneExp, gmt_file, dimension = 1, score_names = emt_name)
      obj@meta.data[[emt_name]] <- nnPCA_Result[rownames(obj@meta.data), emt_name]
      
    } else if (method == "AUCell") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      AUCell_Result <- Execute_AUCell(geneExp, gmt_file, score_names = emt_name)
      obj@meta.data[[emt_name]] <- AUCell_Result[rownames(obj@meta.data), emt_name]
      
    } else if (method == "ssGSEA") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      ssGSVA_Result <- Execute_ssGSVA(geneExp, gmt_file, score_names = emt_name)
      obj@meta.data[[emt_name]] <- ssGSVA_Result[rownames(obj@meta.data), emt_name]
      
    } else if (method == "SCSE") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      SCSE_Result <- Execute_SCSE(geneExp, gmt_file, score_names = emt_name)
      obj@meta.data[[emt_name]] <- SCSE_Result[rownames(obj@meta.data), emt_name]
      
    } else if (method == "JASMINE") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      JAS_Result <- Execute_JAS(geneExp, gmt_file, score_names = emt_name)
      obj@meta.data[[emt_name]] <- JAS_Result[rownames(obj@meta.data), emt_name]
    }
    
    return(obj)
  })
  
  names(seurat_list) <- names(seurat_files)
  return(seurat_list)
}


#' Plot EMT Scores Along Pseudotime from Seurat Objects
#'
#' This function visualizes EMT (epithelial–mesenchymal transition) scores 
#' across pseudotime for multiple Seurat objects. Each Seurat object is treated 
#' as a different condition or sample.
#'
#' @param seurat_list A named list of Seurat objects. Each element should be a Seurat object 
#'   with the EMT score already computed and stored in the metadata.
#' @param col_name Character string. Name of the column in \code{meta.data} containing 
#'  values to plot.
#' @param emt_score_col Character string. Name of the column in \code{meta.data} containing 
#'   EMT scores. Default: \code{"EMT_Score"}.
#'
#' @return A \code{ggplot} object showing smoothed EMT score trajectories along pseudotime 
#'   for each condition.
#'
#' @examples
#' # Assume seurat_list has EMT scores calculated
#' p <- plot_EMT_from_objects(seurat_list, col_name = "Pseudotime", emt_score_col = "EMT_Score")
#' 
#' @export

plot_EMT_from_objects <- function(seurat_list, col_name,
                                  emt_score_col) {
  
  # Merge all Seurat object meta.data
  plot_df <- do.call(rbind, lapply(names(seurat_list), function(name) {
    obj <- seurat_list[[name]]
    df <- obj@meta.data[, c(col_name, emt_score_col), drop = FALSE]
    colnames(df) <- c(col_name, emt_score_col)
    df$Condition <- name
    df <- df[order(df[[col_name]]), ]  
    df
  }))
  
  # Draw smooth curves
  p <- ggplot(plot_df, aes_string(x = col_name, y = emt_score_col, color = "Condition")) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
    theme_classic(base_size = 14) +
    labs(x = col_name, y = emt_score_col, color = "Condition")
  
  return(p)
}
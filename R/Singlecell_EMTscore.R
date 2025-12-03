#' Add EMT Scores to Seurat or SCE Objects
#'
#' This function calculates epithelial mesenchymal transition (EMT) scores 
#' for one or more single-cell objects using a variety of scoring methods. 
#' Supported methods include \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"}, 
#' \code{"ssGSEA"}, \code{"SCSE"}, and \code{"JASMINE"}.
#'
#' @param files A character vector of file paths to Seurat or SCE objects (\code{.rds}
#' or \code{.rda}  files).
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
#' - \strong{nnPCA}: Computes nonlinear PCA based EMT scores.  
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
#' @import Matrix patchwork dplyr slingshot ggplot2 SeuratWrappers DoubletFinder monocle Seurat data.table nsprcomp
#' @export
#' 

add_EMT_score <- function(files, gmt_file = NULL,
                          emt_name, method = c("Seurat","nnPCA","AUCell","ssGSEA","SCSE","JASMINE"),
                          nnPCA_dim = 1) {
  
  method <- match.arg(method)
  
  if (is.null(gmt_file)) stop("gmt_file must be provided.")
  Genesets <- read_gmt(gmt_file)
  
  # add label
  names(files) <- gsub("\\.rds$|\\.rda$", "", basename(files))
  
  obj_list <- lapply(names(files), function(name){
    file <- files[name]
    
    if (grepl("\\.rds$", file)){
      obj <- readRDS(file)
    }
    else if (grepl("\\.rda$", file)){
      e <- new.env()
      load(file, envir = e)
      obj <- e[[ls(e)[1]]]
    }
    else{
      stop("Unsupported file type")
    }
    
    if (inherits(obj, "SingleCellExperiment")){
      message("Converting SCE to Seurat")
      obj <- as.Seurat(obj, data = "logcounts")
    }
    else if (!inherits(obj, "Seurat")){
      stop("Unsupported object type")
    }

    obj <- UpdateSeuratObject(obj)
    emt_genes <- Genesets$gene
    
    if (method == "Seurat") {
      obj <- AddModuleScore(obj, features = list(emt_genes), name = "EMT_Score", ctrl = 5)
      colnames(obj@meta.data)[colnames(obj@meta.data) == "EMT_Score1"] <- emt_name
      
      
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
  
  names(obj_list) <- names(files)
  return(obj_list)
}

#' Add Multiple EMT Scores to Seurat or SCE Objects
#'
#' This function calculates epithelial mesenchymal transition (EMT) related scores 
#' for multiple single-cell objects using multiple gene sets defined in a GMT file. 
#' Supported scoring methods include \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"}, 
#' \code{"ssGSEA"}, \code{"SCSE"}, and \code{"JASMINE"}.
#'
#' @param files A character vector of file paths to Seurat or SCE objects (\code{.rds}
#' or \code{.rda} files).
#' @param gmt_file Character string. Path to the gene set file in GMT format.
#'   The GMT file can contain multiple gene sets; a separate EMT score will be computed for each.
#' @param emt_names Character vector. Names of the EMT score columns to add.
#' @param method Character string. Scoring method to use. Must be one of:
#'   \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"}, \code{"ssGSEA"}, 
#'   \code{"SCSE"}, or \code{"JASMINE"}.
#' @param nnPCA_dim Integer. Dimension to use for nnPCA scoring (default = 1).
#' @param cores Integer. Number of CPU cores to use for parallel computation.
#'
#' @details
#' - \strong{Seurat}: Uses \code{Seurat::AddModuleScore} to compute module scores for each gene set.  
#' - \strong{nnPCA}: Performs nonlinear PCA based scoring for all gene sets in parallel.  
#' - \strong{AUCell}: Calculates AUC-based enrichment scores using parallel execution.  
#' - \strong{ssGSEA}: Computes single-sample GSEA based scores for each gene set.  
#' - \strong{SCSE}: Applies the SCSE algorithm to compute EMT-related activity scores.  
#' - \strong{JASMINE}: Computes EMT gene set scores using the JASMINE algorithm.  
#'
#' For each Seurat object, new columns will be added to the \code{meta.data} slot, 
#' corresponding to each gene set name in the provided GMT file.
#'
#' @return A named list of Seurat objects with updated metadata, each containing
#'   additional columns for all computed EMT-related gene set scores.
#'
#' @import Matrix patchwork dplyr slingshot ggplot2 
#'   SeuratWrappers DoubletFinder monocle Seurat data.table nsprcomp GSA
#' @export
add_EMT_score_multiple <- function(files, gmt_file = NULL,emt_names,
                          method = c("Seurat","nnPCA","AUCell","ssGSEA","SCSE","JASMINE"),
                          nnPCA_dim, cores) {
  
  method <- match.arg(method)
  
  if (is.null(gmt_file)) stop("gmt_file must be provided.")
  Genesets_obj <- GSA.read.gmt(gmt_file)
  
  # add label
  names(files) <- gsub(".rds", "", basename(files))
  
  obj_list <- lapply(names(files), function(name){
    
    file <- files[name]
    if (grepl("\\.rds$", file)) {
      obj <- readRDS(file)
    } else if (grepl("\\.rda$", file)) {
      e <- new.env()
      load(file, envir = e)
      obj <- e[[ls(e)[1]]]
    } else {
      stop("Unsupported file type")
    }
    
    if (inherits(obj, "SingleCellExperiment")) {
      message("Converting SCE → Seurat using as.Seurat() ...")
      obj <- as.Seurat(obj, data = "logcounts")
    } else if (!inherits(obj, "Seurat")) {
      stop("Unsupported object type")
    }
    
    obj <- UpdateSeuratObject(obj)
    feature_lists <- lapply(Genesets_obj$genesets, unlist)
    feature_names <- Genesets_obj$geneset.names
    
    if (method == "Seurat") {
      obj <- AddModuleScore(obj, features = feature_lists, name = emt_names, ctrl = 5)
      for (i in seq_along(emt_names)) {
        old_name <- paste0(emt_names[i], i)
        if (old_name %in% colnames(obj@meta.data)) {
          colnames(obj@meta.data)[colnames(obj@meta.data) == old_name] <- emt_names[i]
        }
      }
      
    } else if (method == "nnPCA") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      nnPCA_Result <- Execute_nnPCA_parallel(geneExp, gmt_file, dimension = 1, cores)
      obj@meta.data <- cbind(obj@meta.data, nnPCA_Result[rownames(obj@meta.data), ])
      
    } else if (method == "AUCell") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      AUCell_Result <- Execute_AUCell_parallel(geneExp, gmt_file, cores)
      obj@meta.data <- cbind(obj@meta.data, AUCell_Result[rownames(obj@meta.data), ])
      
    } else if (method == "ssGSEA") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      ssGSVA_Result <- Execute_ssGSEA_parallel(geneExp, gmt_file, cores)
      obj@meta.data <- cbind(obj@meta.data, ssGSVA_Result[rownames(obj@meta.data), ])
      
    } else if (method == "SCSE") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      SCSE_Result <- Execute_SCSE_parallel(geneExp, gmt_file, cores)
      obj@meta.data <- cbind(obj@meta.data, SCSE_Result[rownames(obj@meta.data), ])
      
    } else if (method == "JASMINE") {
      geneExp <- GetAssayData(obj, assay = "RNA", slot = "data")
      JAS_Result <- Execute_JASMINE_parallel(geneExp, gmt_file, cores)
      obj@meta.data <- cbind(obj@meta.data, JAS_Result[rownames(obj@meta.data), ])
    }
    
    return(obj)
  })
  
  names(obj_list) <- names(files)
  return(obj_list)
}




#' Plot EMT Scores Along Pseudotime from Seurat Objects
#'
#' This function visualizes EMT (epithelial mesenchymal transition) scores 
#' across pseudotime for multiple Seurat objects. Each Seurat object is treated 
#' as a different condition or sample.
#'
#' @param obj_list A named list of Seurat objects. Each element should be a Seurat object 
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
#' # Assume obj_list has EMT scores calculated
#' p <- plot_EMT_from_objects(obj_list, col_name = "Pseudotime", emt_score_col = "EMT_Score")
#' 
#' @export

plot_EMT_from_objects <- function(obj_list, col_name,
                                  emt_score_col) {
  
  # Merge all Seurat object meta.data
  plot_df <- do.call(rbind, lapply(names(obj_list), function(name) {
    obj <- obj_list[[name]]
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

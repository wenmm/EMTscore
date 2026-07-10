#' Add EMT Scores to Seurat or SCE Objects
#'
#' This function calculates epithelial-mesenchymal transition (EMT) scores
#' for one or more single-cell objects using multiple scoring methods.
#' Supported methods include \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"},
#' \code{"GSVA"}, \code{"ssGSEA"}, \code{"SCSE"}, and \code{"JASMINE"}.
#'
#' @param objects A named list of Seurat or SingleCellExperiment objects
#'   that have already been loaded in the R session.
#' @param gmt_file Character string. Path to a GMT file containing EMT gene sets.
#' @param emt_name Character string. Name of the EMT score column to add.
#' @param method Character string. Scoring method to use. One of:
#'   \code{"Seurat"}, \code{"nnPCA"}, \code{"AUCell"},\code{"GSVA"},
#'   \code{"ssGSEA"}, \code{"SCSE"}, or \code{"JASMINE"}.
#' @param nnPCA_dim Integer. Dimension for nnPCA (default = 1).
#'
#' @import ggplot2
#' @return A named list of Seurat objects with a new metadata column containing EMT scores.
#' @examples
#' \donttest{
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' A549_TGFB1 <- eh[["EH10293"]]
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' objects <- list(A549_TGFB1 = A549_TGFB1)
#' res <- add_EMT_score(objects, gmt_file,
#'   emt_name = "EMT_Score",
#'   method = "AUCell"
#' )
#' }
#' @export

add_EMT_score <- function(objects,
                          gmt_file,
                          emt_name = "EMT_Score",
                          method = c("Seurat", "nnPCA", "AUCell", "GSVA", "ssGSEA", "SCSE", "JASMINE"),
                          nnPCA_dim = 1) {
  method <- match.arg(method)
  Genesets <- read_gmt(gmt_file)
  emt_genes <- Genesets$gene

  obj_list <- lapply(names(objects), function(name) {
    obj <- objects[[name]]

    # Convert SCE to Seurat
    if (inherits(obj, "SingleCellExperiment")) {
      message("Converting SCE to Seurat: ", name)
      obj <- as.Seurat(obj, data = "logcounts")
    }
    if (!inherits(obj, "Seurat")) {
      stop("Object ", name, " is not a Seurat or SCE object.")
    }

    obj <- UpdateSeuratObject(obj)

    # get gene expression matrix
    geneExp <- GetAssayData(obj, assay = "RNA", layer = "data")

    if (method == "Seurat") {
      obj <- AddModuleScore(obj, features = list(emt_genes), name = emt_name, ctrl = 5)
      obj[[emt_name]] <- obj[[paste0(emt_name, "1"), drop = TRUE]]
      obj[[paste0(emt_name, "1")]] <- NULL
    } else if (method == "nnPCA") {
      r <- Execute_nnPCA(geneExp, gmt_file, dimension = nnPCA_dim, score_names = emt_name)
      if (nnPCA_dim == 1) {
        obj[[emt_name]] <- r[colnames(obj), emt_name]
      } else if (nnPCA_dim > 1) {
        obj <- AddMetaData(obj, r[colnames(obj), ])
      }
    } else if (method == "AUCell") {
      r <- Execute_AUCell(geneExp, gmt_file, score_names = emt_name)
      obj[[emt_name]] <- r[colnames(obj), emt_name]
    } else if (method == "ssGSEA") {
      r <- Execute_ssGSEA(geneExp, gmt_file, score_names = emt_name)
      obj[[emt_name]] <- r[colnames(obj), emt_name]
    } else if (method == "GSVA") {
      r <- Execute_GSVA(geneExp, gmt_file, score_names = emt_name)
      obj[[emt_name]] <- r[colnames(obj), emt_name]
    } else if (method == "SCSE") {
      r <- Execute_SCSE(geneExp, gmt_file, score_names = emt_name)
      obj[[emt_name]] <- r[colnames(obj), emt_name]
    } else if (method == "JASMINE") {
      r <- Execute_JAS(geneExp, gmt_file, score_names = emt_name)
      obj[[emt_name]] <- r[colnames(obj), emt_name]
    }

    return(obj)
  })

  names(obj_list) <- names(objects)
  return(obj_list)
}


#' Add Multiple EMT Scores to Seurat or SCE Objects
#'
#' Computes EMT-related scores for multiple gene sets defined in a GMT file.
#'
#' @param objects A named list of Seurat or SingleCellExperiment objects already in memory.
#' @param gmt_file Character string. Path to GMT file with multiple gene sets.
#' @param emt_names Character vector. Column names to assign to the resulting EMT scores.
#' @param method Character string. Scoring method.
#' @param nnPCA_dim Integer. Dimension used in nnPCA.
#' @param cores Integer. Number of CPU cores for parallel computation.
#'
#' @import ggplot2
#' @return A named list of Seurat objects with multiple EMT score columns added.
#' @examples
#' \donttest{
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' A549_TNF <- eh[["EH10291"]]
#' A549_EGF <- eh[["EH10292"]]
#' A549_TGFB1 <- eh[["EH10293"]]
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' objects <- list(A549_TGFB1 = A549_TGFB1, A549_EGF = A549_EGF, A549_TNF = A549_TNF)
#' EMscore_result <- add_EMT_score_multiple(objects, gmt_file,
#'   emt_names = c("Escore", "Mscore"),
#'   method = "nnPCA", nnPCA_dim = 1, cores = 2
#' )
#' }
#'
#' @export
add_EMT_score_multiple <- function(objects, gmt_file, emt_names = NULL,
                                   method = c("Seurat", "nnPCA", "AUCell", "ssGSEA", "GSVA", "SCSE", "JASMINE"),
                                   nnPCA_dim = 1, cores = 1) {
  method <- match.arg(method)

  Genesets_obj <- GSA.read.gmt(gmt_file)
  feature_lists <- lapply(Genesets_obj$genesets, unlist)

  if (is.null(emt_names)) {
    message("emt_names not provided, using gene set names from GMT file.")
  }

  obj_list <- lapply(names(objects), function(name) {
    obj <- objects[[name]]

    if (inherits(obj, "SingleCellExperiment")) {
      message("Converting SCE to Seurat: ", name)
      obj <- as.Seurat(obj, data = "logcounts")
    }
    if (!inherits(obj, "Seurat")) {
      stop("Object ", name, " is not a Seurat or SCE object.")
    }

    obj <- UpdateSeuratObject(obj)
    geneExp <- GetAssayData(obj, assay = "RNA", layer = "data")

    if (method == "Seurat") {
      obj <- AddModuleScore(obj, features = feature_lists, name = emt_names, ctrl = 5)

      # Rename Seurat default column names
      for (i in seq_along(emt_names)) {
        old <- paste0(emt_names[i], i)
        if (old %in% colnames(obj[[]])) {
          obj[[emt_names[i]]] <- obj[[old, drop = TRUE]]
          obj[[old]] <- NULL
        }
      }
    } else if (method == "nnPCA") {
      r <- Execute_nnPCA_parallel(geneExp, gmt_file, dimension = nnPCA_dim, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    } else if (method == "AUCell") {
      r <- Execute_AUCell_parallel(geneExp, gmt_file, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    } else if (method == "ssGSEA") {
      r <- Execute_ssGSEA_parallel(geneExp, gmt_file, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    } else if (method == "GSVA") {
      r <- Execute_GSVA_parallel(geneExp, gmt_file, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    } else if (method == "SCSE") {
      r <- Execute_SCSE_parallel(geneExp, gmt_file, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    } else if (method == "JASMINE") {
      r <- Execute_JASMINE_parallel(geneExp, gmt_file, cores = cores)
      if (!is.null(emt_names)) {
        colnames(r) <- emt_names
      }
      obj <- AddMetaData(obj, as.data.frame(as.matrix(r)))
    }

    return(obj)
  })

  names(obj_list) <- names(objects)
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
#' \donttest{
#' library(ExperimentHub)
#' eh <- ExperimentHub()
#' A549_TGFB1 <- eh[["EH10293"]]
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' objs <- add_EMT_score(list(A549_TGFB1 = A549_TGFB1), gmt_file,
#'   emt_name = "EMT_Score", method = "AUCell"
#' )
#' p <- plot_EMT_from_objects(objs,
#'   col_name = "nCount_RNA",
#'   emt_score_col = "EMT_Score"
#' )
#' }
#' @export

plot_EMT_from_objects <- function(obj_list, col_name,
                                  emt_score_col) {
  # Merge all Seurat object meta.data
  plot_df <- do.call(rbind, lapply(names(obj_list), function(name) {
    obj <- obj_list[[name]]
    df <- obj[[c(col_name, emt_score_col)]]
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

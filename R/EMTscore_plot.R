#' Prepare data for plotting
#'
#' Merge cell annotation and score result, and set rownames based on a user-specified column.
#'
#' @param cell_annotation_file A data.frame containing cell annotation.
#' @param score_result A data.frame or matrix of computed scores (nnPCA, AUCell, etc.).
#' @param merge_colname Character, the column of cell_annotation_file to use as row names.
#'
#' @return A data.frame ready for plotting.
#' @examples
#' data(cell_annotation_file)
#' data(nnPCA_Result_multiple)
#' data_for_plot <- data_prepare(cell_annotation_file, nnPCA_Result_multiple, merge_colname = "name")
#' @export
data_prepare <- function(cell_annotation_file, score_result, merge_colname) {
  merge.all <- function(x, ..., by = "row.names") {
    L <- list(...)
    for (i in seq_along(L)) {
      x <- merge(x, L[[i]], by = by)
      rownames(x) <- x$Row.names
      x$Row.names <- NULL
    }
    return(x)
  }
  rownames(cell_annotation_file) <- cell_annotation_file[[merge_colname]]
  data_for_analysis <- merge.all(cell_annotation_file, score_result)
  return(data_for_analysis)
}

#' Plot E vs M scores
#'
#' Scatter plot of E vs M scores with mean +/- SD per cell type.
#'
#' @param data_for_plot A data.frame prepared by plot_data_prepare.
#' @param E_colname Character, column name for E score.
#' @param M_colname Character, column name for M score.
#' @param celltype_colname Character, column name for cell type annotation.
#' @param colors Vector of colors for cell types.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @examples
#' data(data_for_plot)
#' plot <- Execute_E_M_plot(
#'   data_for_plot,
#'   E_colname = "Panchy_et_al_E_signature",
#'   M_colname = "Panchy_et_al_M_signature",
#'   celltype_colname = "celltype_annotation",
#'   colors = c(
#'     "#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
#'     "#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6"
#'   )
#' )
#' @export
Execute_E_M_plot <- function(data_for_plot,
                             E_colname,
                             M_colname,
                             celltype_colname,
                             colors = c(
                               "#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
                               "#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6"
                             )) {
  E <- sym(E_colname)
  M <- sym(M_colname)
  CellType <- sym(celltype_colname)

  data_for_plot[[E_colname]] <- as.numeric(unlist(data_for_plot[[E_colname]]))
  data_for_plot[[M_colname]] <- as.numeric(unlist(data_for_plot[[M_colname]]))


  data_for_plot <- data_for_plot %>%
    dplyr::filter(
      !is.na(!!E),
      !is.na(!!M)
    )

  p <- ggplot(data_for_plot, aes(x = !!E, y = !!M, fill = !!CellType)) +
    geom_point(size = 2.5, alpha = 0.3, aes(color = !!CellType)) +
    scale_colour_manual(values = colors) +
    stat_density2d(aes(x = !!E, y = !!M), bins = 10, alpha = 0.2, geom = "polygon")


  sbg <- data_for_plot %>%
    group_by(!!CellType) %>%
    summarise(
      count = n(),
      mE = mean(as.numeric(unlist(!!E)), na.rm = TRUE),
      sdE = sd(as.numeric(unlist(!!E)), na.rm = TRUE),
      mM = mean(as.numeric(unlist(!!M)), na.rm = TRUE),
      sdM = sd(as.numeric(unlist(!!M)), na.rm = TRUE),
      .groups = "drop"
    )


  sbg <- sbg[complete.cases(sbg), ]

  p2 <- p +
    geom_errorbar(
      data = sbg,
      aes(
        x = mE, y = mM,
        ymin = mM - 0.8 * sdM,
        ymax = mM + 0.8 * sdM
      ),
      width = 0, size = 2
    ) +
    geom_errorbarh(
      data = sbg,
      aes(
        x = mE, y = mM,
        xmin = mE - 0.8 * sdE,
        xmax = mE + 0.8 * sdE
      ),
      height = 0, size = 2
    ) +
    geom_point(
      data = sbg,
      aes(x = mE, y = mM, fill = !!CellType),
      color = "white", shape = 21, size = 6,
      alpha = 1, show.legend = FALSE, stroke = 1.5
    ) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    theme(
      panel.border = element_rect(color = "black", size = 1.5, fill = NA),
      title = element_text(
        size = 12, color = "black",
        face = "italic", hjust = 0.5, lineheight = 0.2
      ),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks.length.x = unit(0.15, "cm"),
      axis.ticks.length.y = unit(0.15, "cm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.box.background = element_rect(fill = NA, color = "black", linetype = 1)
    )

  return(p2)
}

#' Plot M dimension scores
#'
#' Scatter plot of two M scores with mean +/- SD per cell type.
#'
#' @param data_for_plot A data.frame prepared by plot_data_prepare.
#' @param M1_colname Column name for M1 score.
#' @param M2_colname Column name for M2 score.
#' @param celltype_colname Column name for cell type annotation.
#' @param colors Vector of colors for cell types.
#' @import ggplot2
#' @return A ggplot object.
#' @examples
#' data(cell_annotation_file)
#' data(geneExp)
#' gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
#' nnPCA_Mscore <- Execute_nnPCA(geneExp, gmt_file,
#'   dimension = 2, score_names = c("M1_score", "M2_score"))
#' data_for_plot <- data_prepare(cell_annotation_file, nnPCA_Mscore, merge_colname = "name")
#' plot2 <- Execute_M_dimension_plot(
#'   data_for_plot,
#'   M1_colname = "M1_score",
#'   M2_colname = "M2_score",
#'   celltype_colname = "celltype_annotation",
#'   colors = c(
#'     "#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
#'     "#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6"
#'   )
#' )
#'
#' @export
Execute_M_dimension_plot <- function(data_for_plot,
                                     M1_colname,
                                     M2_colname,
                                     celltype_colname,
                                     colors = c(
                                       "#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
                                       "#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6"
                                     )) {
  M1 <- sym(M1_colname)
  M2 <- sym(M2_colname)
  CellType <- sym(celltype_colname)

  p <- ggplot(data_for_plot, aes(x = !!M1, y = !!M2, fill = !!CellType)) +
    geom_point(size = 2.5, alpha = 0.3, aes(color = !!CellType)) +
    scale_colour_manual(values = colors) +
    stat_density2d(aes(x = !!M1, y = !!M2), bins = 10, alpha = 0.2, geom = "polygon")

  sbg <- data_for_plot %>%
    group_by(!!CellType) %>%
    summarise(
      count = n(),
      mM1 = mean(!!M1, na.rm = TRUE),
      sdM1 = sd(!!M1, na.rm = TRUE),
      mM2 = mean(!!M2, na.rm = TRUE),
      sdM2 = sd(!!M2, na.rm = TRUE)
    )

  p2 <- p +
    geom_errorbar(
      data = sbg,
      aes(
        x = mM1, y = mM2,
        ymin = mM2 - 0.8 * sdM2,
        ymax = mM2 + 0.8 * sdM2
      ),
      width = 0, size = 2
    ) +
    geom_errorbarh(
      data = sbg,
      aes(
        x = mM1, y = mM2,
        xmin = mM1 - 0.8 * sdM1,
        xmax = mM1 + 0.8 * sdM1
      ),
      height = 0, size = 2
    ) +
    geom_point(
      data = sbg,
      aes(x = mM1, y = mM2, fill = !!CellType),
      color = "white", shape = 21, size = 6,
      alpha = 1, show.legend = FALSE, stroke = 1.5
    ) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    theme(
      panel.border = element_rect(color = "black", size = 1.5, fill = NA),
      title = element_text(
        size = 12, color = "black",
        face = "italic", hjust = 0.5, lineheight = 0.2
      ),
      axis.text = element_text(size = 12, color = "black"),
      axis.ticks.length.x = unit(0.15, "cm"),
      axis.ticks.length.y = unit(0.15, "cm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.box.background = element_rect(fill = NA, color = "black", linetype = 1)
    )

  return(p2)
}

#' Arrange multiple ggplot objects with subtitles and common legend
#'
#' @param plots_list List of ggplot objects.
#' @param ncol_per_row Number of plots per row.
#' @param subtitles Optional vector of subtitles (bold).
#' @param fig_title Optional main figure title.
#' @param common_legend Logical, whether to use a common legend.
#' @param legend_position Position of the legend ("bottom", "right", etc.).
#'
#' @return Combined ggplot object.
#' @importFrom ggpubr ggarrange annotate_figure text_grob
#' @export
Arrange_plots <- function(plots_list,
                          ncol_per_row = 2,
                          subtitles = NULL,
                          fig_title = NULL,
                          common_legend = TRUE,
                          legend_position = "bottom") {
  if (!is.null(subtitles)) {
    if (length(subtitles) != length(plots_list)) stop("Length of 'subtitles' must match number of plots")
    for (i in seq_along(plots_list)) {
      plots_list[[i]] <- plots_list[[i]] +
        labs(subtitle = subtitles[i]) +
        theme(plot.subtitle = element_text(face = "bold", hjust = 0.5))
    }
  }

  combined <- ggarrange(
    plotlist = plots_list,
    ncol = ncol_per_row,
    nrow = ceiling(length(plots_list) / ncol_per_row),
    common.legend = common_legend,
    legend = legend_position
  )

  if (!is.null(fig_title)) {
    combined <- annotate_figure(combined,
      top = text_grob(fig_title, face = "bold")
    )
  }

  return(combined)
}

############################# Execution of Plot function

## Function Input:  gene expression matrix and cell type information

## Function Output: Figure of heat map
#' present result of heat map of genes that contribute to M score of dimension 1 and disension 2
#'
#' @param geneExp gene expression matrix
#' @param geneList_M M signature gene list in dataframe, we supply example list such as geneList_M, M_signature_for_cancer, M_signature_for_cell
#'
#' @return Figure represent EMT score for different methods and gene sets
#' @examples
#' data(geneExp)
#' data(Tan_et_al_cell_line_M_signature)
#' geneList_M <- Tan_et_al_cell_line_M_signature
#' p <- plot_heatmap_function(t(geneExp), geneList_M)
#' @export
#'

plot_heatmap_function <- function(geneExp, geneList_M) {
  geneExp_M <- geneExp[, colnames(geneExp) %in% geneList_M$GeneName]
  pc_feature_M <- nsprcomp(as.matrix(geneExp_M), nneg = TRUE, ncomp = 2)

  # Get the top genes based on pcs
  result <- data.frame(pc_feature_M$rotation)
  PC1_gname <- rownames(result[order(-result$PC1), ])[seq_len(10)]
  PC2_gname <- rownames(result[order(-result$PC2), ])[seq_len(10)]

  # Extract expression data for top genes
  PC1_gname_Exp <- geneExp[, colnames(geneExp) %in% PC1_gname]
  PC2_gname_Exp <- geneExp[, colnames(geneExp) %in% PC2_gname]

  # Arrange data for heatmap
  new_df1 <- t(PC1_gname_Exp)[PC1_gname, ]
  new_df2 <- t(PC2_gname_Exp)[PC2_gname, ]

  result_plot <- rbind(new_df1, new_df2)

  # Generate M scores for each component
  M_score <- data.frame(pc_feature_M$x)
  colnames(M_score) <- c("M_PC1_score", "M_PC2_score")

  # Reorder expression data to match M scores
  sorted_data_M <- M_score[colnames(result_plot), ]


  # PCs and annotations
  PC_gname <- append(PC1_gname, PC2_gname)
  result_plot2 <- result[PC_gname, ]

  # Define color palettes: white as the minimum and the darker color as the max
  expr_colors <- colorRamp2(
    c(0, 5, 15),
    c("#2166AC", "white", "#c21f30")
  )

  pc1_max <- max(result_plot2$PC1, na.rm = TRUE)

  pc1_colors <- colorRamp2(
    c(0, pc1_max),
    c("white", "#2166AC")
  )

  pc2_max <- max(result_plot2$PC2, na.rm = TRUE)
  pc2_colors <- colorRamp2(
    c(0, pc2_max),
    c("white", "#c21f30")
  )

  m_pc1_colors <- colorRamp2(
    c(min(sorted_data_M$M_PC1_score), max(sorted_data_M$M_PC1_score)),
    c("white", "#322f3b")
  )
  m_pc2_colors <- colorRamp2(
    c(min(sorted_data_M$M_PC2_score), max(sorted_data_M$M_PC2_score)),
    c("white", "#9b1e64")
  )

  # Categorical palette
  label_colors <- c("#2166AC", "#c21f30")
  names(label_colors) <- c("M_PC1", "M_PC2")

  # Define annotations
  row_ha <- rowAnnotation(
    PC1 = result_plot2$PC1, PC2 = result_plot2$PC2,
    col = list(
      PC1 = pc1_colors,
      PC2 = pc2_colors
    )
  )


  # Define categorical color
  row_labels <- c(rep("M_PC1", length(PC1_gname)), rep("M_PC2", length(PC2_gname)))
  row_label_annotation <- rowAnnotation(
    Label = row_labels,
    col = list(Label = label_colors)
  )

  # Define continuous color map for M scores
  ha <- HeatmapAnnotation(
    M_PC1_score = sorted_data_M$M_PC1_score,
    M_PC2_score = sorted_data_M$M_PC2_score,
    col = list(
      M_PC1_score = m_pc1_colors,
      M_PC2_score = m_pc2_colors
    )
  )

  # Generate heatmap with specified color schemes
  p <- Heatmap(
    result_plot,
    col = expr_colors,
    cluster_rows = FALSE,
    top_annotation = ha,
    right_annotation = row_ha,
    left_annotation = row_label_annotation,
    show_column_names = FALSE,
    name = "Gene Expr",
  )
  return(p)
}

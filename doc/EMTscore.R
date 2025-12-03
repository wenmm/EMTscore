## ----setup, include=FALSE-----------------------------------------------------
# Load required libraries
library(nsprcomp)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(curl)
library(ggpubr)
library(pheatmap)
library(grid)
library(circlize)
library(paletteer)
library(ggthemes)
library(ComplexHeatmap)
library(EMTscore)
library(AUCell)
library(DoubletFinder)
library(SeuratWrappers)
library(monocle)
library(slingshot)
library(Seurat)
library(GSA)
library(mclust)
library(ggalluvial)
library(RColorBrewer)
library(hrbrthemes)
library(Cairo)
options(mc.cores = 4)
options(bitmapType='cairo')


# Optional: set global chunk options
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)

## -----------------------------------------------------------------------------
data(cell_annotation_file)
head(cell_annotation_file)

## -----------------------------------------------------------------------------
url <- "https://zenodo.org/records/17353682/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
head(geneExp[, 1:5], 4)
dim(geneExp)

## -----------------------------------------------------------------------------
data("Panchy_et_al_E_signature", package = "EMTscore")
data("Panchy_et_al_M_signature", package = "EMTscore")
head(Panchy_et_al_M_signature)

## -----------------------------------------------------------------------------
gene_sets <- list(
Panchy_et_al_E_signature = Panchy_et_al_E_signature$GeneName,
Panchy_et_al_M_signature = Panchy_et_al_M_signature$GeneName
)
write_gmt(gene_sets, "EM_signature.gmt")

## -----------------------------------------------------------------------------
nnPCA_Result_multiple <- Execute_nnPCA_parallel(geneExp, 'EM_signature.gmt', dimension = 1, cores = 2)
AUCell_Result_multiple <- Execute_AUCell_parallel(geneExp, 'EM_signature.gmt', cores = 2)
ssGSEA_Result_multiple <- Execute_ssGSEA_parallel(geneExp, 'EM_signature.gmt', cores = 2)
JASMINE_Result_multiple <- Execute_JASMINE_parallel(geneExp, 'EM_signature.gmt', cores = 2)
SCSE_Result_multiple <- Execute_SCSE_parallel(geneExp, 'EM_signature.gmt', cores = 2)

## -----------------------------------------------------------------------------
data_for_plot <- data_prepare(cell_annotation_file, nnPCA_Result_multiple, merge_colname = "name")

## -----------------------------------------------------------------------------
plot1 <- Execute_E_M_plot(
data_for_plot,
E_colname = "Panchy_et_al_E_signature",
M_colname = "Panchy_et_al_M_signature",
celltype_colname = "celltype_annotation",
colors = c("#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
"#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6")
)
plot1

## -----------------------------------------------------------------------------
# Prepare M gene set and compute nnPCA

gene_set <- list(Panchy_et_al_M_signature = Panchy_et_al_M_signature$GeneName)
write_gmt(gene_set, "M_signature.gmt")

nnPCA_Mscore <- Execute_nnPCA(geneExp, "M_signature.gmt", dimension=2, score_names=c('M1_score','M2_score'))
data_for_plot <- data_prepare(cell_annotation_file, nnPCA_Mscore, merge_colname = "name")

plot2 <- Execute_M_dimension_plot(
data_for_plot,
M1_colname = "M1_score",
M2_colname = "M2_score",
celltype_colname = "celltype_annotation",
colors = c("#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
"#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6")
)
plot2

## -----------------------------------------------------------------------------
combined_plot <- Arrange_plots(
plots_list = list(plot1, plot2),
ncol_per_row = 2,
subtitles = c("E vs M", "M1 vs M2"),
fig_title = "Panchy_et_al"
)
combined_plot

## -----------------------------------------------------------------------------
colors = c("#F87189", "#CE9031", "#A48CF5", "#97A430", "#39A7D0", "#E57D5F",
"#84C7B9", "#E1AF64", "#C26CCF", "#B0BF43", "#57C3E8", "#F29D9E", "#92AAE6")
p_hist <- data_for_plot %>%
     ggplot( aes(x=M1_score, fill=celltype_annotation)) +
     geom_histogram(alpha=0.6, position = 'identity') +
     scale_fill_manual(values=colors) +
     theme_ipsum() +
     labs(fill="")
p_hist

## -----------------------------------------------------------------------------
# Example heatmap

plot_heatmap_function(t(geneExp), Panchy_et_al_M_signature)


## -----------------------------------------------------------------------------
obj <- system.file("extdata", "example_seurat_obj.rds", package = "EMTscore")
files <- c("/home/hwen6/project/EMT/data/GSE147405Cook/A549_TGFB1.rds","/home/hwen6/project/EMT/data/GSE147405Cook/A549_TNF.rds", "/home/hwen6/project/EMT/data/GSE147405Cook/A549_EGF.rds")
#files <- c(obj)
gmt <- system.file("extdata", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.gmt", package = "EMTscore")

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "AUCell",nnPCA_dim = 1)
p_AUCell <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_AUCell

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "SCSE")
p_SCSE <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_SCSE

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "ssGSEA")
p_ssGSEA <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_ssGSEA

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "Seurat")
p_Seurat <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_Seurat

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "nnPCA")
p_nnPCA <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_nnPCA

seurat_objs <- add_EMT_score(files, gmt_file = gmt, emt_name = "EMT_score", method = "JASMINE")
p_JASMINE <- plot_EMT_from_objects(seurat_objs, col_name = "Pseudotime", emt_score_col = "EMT_score")
p_JASMINE

## ----message=FALSE------------------------------------------------------------
gmt_file <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")
files <- c("/home/hwen6/project/EMT/data/GSE147405Cook/A549_TGFB1.rds","/home/hwen6/project/EMT/data/GSE147405Cook/A549_TNF.rds", "/home/hwen6/project/EMT/data/GSE147405Cook/A549_EGF.rds")

emt_name = c("Escore", "Mscore")
result <- add_EMT_score_multiple(files, gmt_file, emt_name, method = "Seurat", nnPCA_dim = 1, cores = 10)

#save into pdf format

plot_all_clusters <- function(result, method = c("Kmeans", "GMM"),emt_name,  n_clusters = 3) {
  method <- match.arg(method)
  
  for (name in names(result)) {
    cat("Processing:", name, "\n")
    obj <- result[[name]]
    
    # Extract E and M signature data
    sig_df <- obj@meta.data[, emt_name, drop = FALSE]
    colnames(sig_df) <- c("Escore", "Mscore")
    
    # Run clustering
    if (method == "GMM") {
      cl <- predict_cluster_labels(sig_df, method = "GMM", PC_name = c("Escore", "Mscore"))
    } else if (method == "Kmeans") {
      cl <- predict_cluster_labels(sig_df, method = "Kmeans", n_clusters = n_clusters, PC_name = c("Escore", "Mscore"))
    }
    
    # Scatter plot data
    scatter_df <- data.frame(sig_df, Cluster = as.factor(cl))
    p1 <- ggplot(scatter_df, aes(x = Escore, y = Mscore, color = Cluster)) +
      geom_point(size = 1.5, alpha = 0.8) +
      scale_color_brewer(palette = "Set1") +
      theme_classic(base_size = 14) +
      labs(title = paste0(name, " - ", method, " clustering"), x = "Escore", y = "Mscore")
    
    # Save scatter plot
    scatter_file <- paste0(name, "_EM_", tolower(method), ".pdf")
    ggsave(scatter_file, plot = p1, width = 6, height = 5)
    message("Sved scatter plot: ", scatter_file)
    
    # Sankey data
    cl_df <- data.frame(cell = rownames(sig_df), Cluster = as.character(unlist(cl)))
    true_df <- data.frame(cell = rownames(obj@meta.data),
                          TrueLabel = as.character(unlist(obj@meta.data$Time)))
    df_merge <- inner_join(cl_df, true_df, by = "cell")
    df_count <- df_merge %>%
      group_by(Cluster, TrueLabel) %>%
      summarise(Freq = n(), .groups = "drop")
    
    # Sankey plot
    p2 <- ggplot(df_count, aes(axis1 = Cluster, axis2 = TrueLabel, y = Freq)) +
      geom_alluvium(aes(fill = Cluster), width = 1/12) +
      geom_stratum(width = 1/8, fill = "grey90", color = "black") +
      geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c(method, "TrueLabel"), expand = c(.1, .1)) +
      scale_fill_brewer(palette = "Set1") +
      theme_minimal(base_size = 14) +
      labs(title = paste("Sankey Diagram:", method, "Cluster vs True Cell Labels -", name),
           y = "Number of Cells", x = "")
    
    # Save Sankey plot
    sankey_file <- paste0(name, "_EM_Sankey_Diagram_", tolower(method), ".pdf")
    ggsave(sankey_file, plot = p2, width = 7, height = 5)
    message("Saved Sankey plot: ", sankey_file)
  }
  
  message("All plots successfully generated and saved.")
}

plot_all_clusters(result, method = "Kmeans", emt_name, n_clusters = 3)
# or
plot_all_clusters(result, method = "GMM", emt_name)


## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
#SCSE_Result_multiple <- Execute_SCSE_parallel(geneExp,'/home/hwen6/project/EMTscore/2025_10_29/EMTscore/inst/extdata/c2.all.v2025.1.Hs.symbols.gmt', cores = 20)
data(SCSE_Result_EMT)
data(SCSE_Result_multiple)
result <- correlate_sample_scores(score_mat1 = SCSE_Result_multiple, score_mat2 = SCSE_Result_EMT, method = "spearman")
head(result)

## -----------------------------------------------------------------------------
signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
scores <- compute_Signature_score(geneExp, signature_file , score_name = "stemness_score")
head(scores)

## -----------------------------------------------------------------------------
signature_file <- system.file("extdata", "cellular_senescence_sig.tsv", package = "EMTscore")
scores <- compute_Signature_score(geneExp, signature_file, score_name = "senescence_score")
head(scores)

## -----------------------------------------------------------------------------
files <- c("/home/hwen6/project/EMT/data/GSE147405Cook/A549_TGFB1.rds","/home/hwen6/project/EMT/data/GSE147405Cook/A549_TNF.rds", "/home/hwen6/project/EMT/data/GSE147405Cook/A549_EGF.rds")

#obj <- system.file("extdata", "example_seurat_obj.rds", package = "EMTscore")
#files <- c(obj)
signature_file1 <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
signature_file2 <- system.file("extdata", "cellular_senescence_sig.tsv", package = "EMTscore")

signature_files <- c(signature_file1, signature_file2)

result <- compute_Signature_score_SingleCell(files, signature_files, score_name = c("Stemness_Score", "Senescence_Score"))
head(result$A549_TGFB1@meta.data)

## -----------------------------------------------------------------------------

gmt_file = "EM_signature.gmt"

files <- c("/home/hwen6/project/EMT/data/GSE147405Cook/A549_TGFB1.rds","/home/hwen6/project/EMT/data/GSE147405Cook/A549_TNF.rds", "/home/hwen6/project/EMT/data/GSE147405Cook/A549_EGF.rds")

EMscore_result <- add_EMT_score_multiple(files, gmt_file,emt_name = c("Escore", "Mscore"),
                          method = "nnPCA",
                          nnPCA_dim = 1, cores = 2)

head(EMscore_result$A549_TGFB1@meta.data)

## -----------------------------------------------------------------------------
df <- result$A549_TGFB1@meta.data
cor_test <- cor.test(df$Senescence_Score, df$Stemness_Score, method = "pearson")
R_val <- round(cor_test$estimate, 2)
p_val <- format(cor_test$p.value, scientific = TRUE, digits = 2)
n_clust <- length(unique(df$Time))
palette_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_clust)
ggplot(df, aes(x = Senescence_Score, y = Stemness_Score, color = Time)) +
     geom_point(size = 2.5, alpha = 0.8) +
     geom_smooth(method = "lm", se = FALSE, size = 1) +
     scale_color_manual(values = palette_colors) +
     annotate("text",
              x = min(df$Senescence_Score),
              y = max(df$Stemness_Score),
              label = paste0("R = ", R_val, ", p < ", p_val),
              hjust = 0, size = 5) +
    theme_bw() +
  theme(
    panel.grid = element_blank(),      
    panel.border = element_rect(       
      colour = "black", 
      fill = NA, 
      size = 1
    )) + theme(
  text = element_text(colour = "black")
) + theme(
  axis.title = element_text(size = 14),     
  axis.text  = element_text(size = 14),     
  legend.title = element_text(size = 13),   
  legend.text = element_text(size = 11),    
  plot.title = element_text(size = 12, face = "bold") 
)
# 5.25 3.84

## -----------------------------------------------------------------------------
# Extract metadata
df1 <- result$A549_TGFB1@meta.data
df2 <- EMscore_result$A549_TGFB1@meta.data

# Add cell names for merging
df1$cell <- rownames(df1)
df2$cell <- rownames(df2)

# Merge metadata
df <- df1 %>%
  select(cell, Stemness_Score, Senescence_Score, Time) %>%
  left_join(df2 %>% select(cell, 
                           Panchy_et_al_E_signature, 
                           Panchy_et_al_M_signature), by = "cell")

# Automatically generate color palette for Time groups
n_groups <- length(unique(df$Time))
palette_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_groups)

# Define combinations to plot
x_vars <- c("Stemness_Score", "Senescence_Score")
y_vars <- c("Panchy_et_al_E_signature", "Panchy_et_al_M_signature")

# Loop through all combinations
for (x_var in x_vars) {
  for (y_var in y_vars) {
    
    # Pearson correlation
    cor_test <- cor.test(df[[x_var]], df[[y_var]])
    R_val <- round(cor_test$estimate, 2)
    p_val <- format(cor_test$p.value, scientific = TRUE, digits = 2)
    
    # Plot
    p <- ggplot(df, aes_string(x = x_var, y = y_var, color = "Time")) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = FALSE, size = 1) +
      scale_color_manual(values = palette_colors) +
      annotate("text",
               x = min(df[[x_var]]),
               y = max(df[[y_var]]),
               label = paste0("R = ", R_val, ", p < ", p_val),
               hjust = 0,
               size = 5) +
      labs(
        x = gsub("_", " ", x_var),
        y = gsub("_", " ", y_var),
        title = paste(gsub("_", " ", y_var), "vs", gsub("_", " ", x_var))
      ) +
      theme_bw() +
  theme(
    panel.grid = element_blank(),      
    panel.border = element_rect(      
      colour = "black", 
      fill = NA, 
      size = 1
    )) + theme(
  text = element_text(colour = "black")
) + theme(
  axis.title = element_text(size = 14),     
  axis.text  = element_text(size = 14),    
  legend.title = element_text(size = 13),   
  legend.text = element_text(size = 11),    
  plot.title = element_text(size = 12, face = "bold")  
)
    
    print(p)
  }
}

## -----------------------------------------------------------------------------
sessionInfo()


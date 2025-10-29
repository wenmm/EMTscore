
# Generate synthetic count data representing a single lineage (1000 cells)

set.seed(123)  # for reproducibility

means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})


load("/home/hwen6/project/EMTscore/2025_09_24/EMTscore/vignettes/merged_genes.rda")
rownames(counts) <- merged_genes
colnames(counts) <- paste0('c',1:300)

library(scater)

sce <- SingleCellExperiment(assays = list(counts = counts))

geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
library(Seurat)

# 确保 assays 存在
sce <- logNormCounts(sce)
if (!"logcounts" %in% assayNames(sce)) {
  logcounts(sce) <- log1p(counts(sce))
}

# 提取矩阵
counts_mat <- as.matrix(assay(sce, "counts"))
logcounts_mat <- as.matrix(assay(sce, "logcounts"))

# --- 构建 Seurat v5 RNA assay 对象 ---
seurat_obj <- CreateSeuratObject(counts = counts_mat, assay = "RNA")

seurat_obj <- SetAssayData(
  object   = seurat_obj,
  assay    = "RNA",
  layer    = "data",      # equivalent to normalized data
  new.data = logcounts_mat
)

seurat_obj <- NormalizeData(seurat_obj)

# 检查结果
print(Assays(seurat_obj))
print(Layers(seurat_obj[["RNA"]]))


saveRDS(seurat_obj, file = "example_seurat_obj.rds")
saveRDS(sce, file = "example_sce.rds")

colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == "slingPseudotime_1"] <- "Pseudotime"


files <- c("/home/hwen6/project/single_cell/collaborate2/Kwende_et_al_manuscript_data_20241114/Kwende_et_al_manuscript_data/joshi_seuratobj.rds", "/home/hwen6/project/single_cell/collaborate2/Kwende_et_al_manuscript_data_20241114/Kwende_et_al_manuscript_data/joshi_seuratobj.rds")
files <- c("/home/hwen6/project/EMTscore/2025_09_24/EMTscore/vignettes/example_seurat_obj.rds", "/home/hwen6/project/EMTscore/2025_09_24/EMTscore/vignettes/example_seurat_obj.rds")

gmt <- "/home/hwen6/project/EMTscore/multiple_genesets/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.gmt"
gmt <- "/home/hwen6/project/EMTscore/2025_09_24/EMTscore/vignettes/mouse_emt.gmt"
seurat_objs <- add_EMT_score(files, gmt_file = gmt, method = "JASMINE")

plot_EMT_from_objects(seurat_objs,col = "Pseudotime",
                      emt_score_col = "EMT_Score")




library(Seurat)

# Suppose geneExp is your normalized matrix (genes x cells)
# Make sure rownames = genes, colnames = cells
geneExp <- as.matrix(geneExp)

# Create a minimal Seurat object (dummy counts just for structure)
seurat_obj <- CreateSeuratObject(
  counts = geneExp * 0 + 1,  # placeholder matrix
  project = "SCLC",
  min.cells = 0,
  min.features = 0
)

# Replace the normalized data slot properly (Seurat v5 uses SetAssayData)
seurat_obj <- SetAssayData(
  object = seurat_obj,
  slot = "data",
  new.data = geneExp,
  assay = "RNA"
)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj


sce <- SingleCellExperiment(assays = list(counts = geneExp))


FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

library(mclust, quietly = TRUE)
cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

cl2 <- kmeans(rd1, centers = 4)$cluster
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

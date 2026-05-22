#' Predict Cluster Labels from PCA Results Using GMM or Kmeans
#'
#' This function clusters PCA results using Gaussian Mixture Model (GMM) or K-means,
#' identifies the clusters with maximum Escore and Mscore, assigns labels "E", "M",
#' and "EM1", "EM2", ... to other clusters, and returns a named vector of final cluster
#' labels for each sample.
#'
#' @param pca_result A numeric matrix or data frame of PCA results (rows = samples, columns = PCs).
#' @param method Character, either "GMM" (Gaussian Mixture Model) or "Kmeans" for clustering. Default is "GMM".
#' @param n_clusters Numeric, number of clusters for Kmeans. Ignored if method = "GMM".
#' @param PC_name Character vector of length 2 specifying the column names for Escore and Mscore. Default is c("Escore", "Mscore").
#' @import mclust
#'
#' @return A named character vector of cluster labels for each sample.
#' @examples
#' data(nnPCA_Result_multiple)
#' labels_gmm <- predict_cluster_labels(nnPCA_Result_multiple, method = "GMM", n_clusters = 3, PC_name = c("Escore", "Mscore"))
#' labels_km  <- predict_cluster_labels(nnPCA_Result_multiple, method = "Kmeans", n_clusters = 3,PC_name = c("Escore", "Mscore"))
#' @export
predict_cluster_labels <- function(pca_result, method = "GMM", n_clusters = 3, PC_name = c("Escore", "Mscore")) {
  
  if (!(method %in% c("GMM", "Kmeans"))) stop("method must be 'GMM' or 'Kmeans'")
  
  if (method == "GMM") {
    cluster_result <- Mclust(pca_result, G = n_clusters)
    cl1 <- cluster_result$classification
    cluster_center_Mean <- cluster_result$parameters$mean
    rownames(cluster_center_Mean) <- PC_name
  } else { # Kmeans
    cluster_result <- kmeans(pca_result, centers = n_clusters)
    cl1 <- cluster_result$cluster
    cluster_center_Mean <- cluster_result$centers
    colnames(cluster_center_Mean) <- PC_name
  }
  
  # Determine E and M clusters
  if (method == "GMM") {
    E_cluster <- which.max(cluster_center_Mean["Escore", ])
    M_cluster <- which.max(cluster_center_Mean["Mscore", ])
    all_clusters <- unique(cl1)
  } else {
    E_cluster <- which.max(cluster_center_Mean[, "Escore"])
    M_cluster <- which.max(cluster_center_Mean[, "Mscore"])
    all_clusters <- rownames(cluster_center_Mean)
    if (is.null(all_clusters)) all_clusters <- as.character(1:nrow(cluster_center_Mean))
  }
  
  other_clusters <- setdiff(all_clusters, c(E_cluster, M_cluster))
  
  # Assign labels
  cluster_labels <- character(length(all_clusters))
  names(cluster_labels) <- all_clusters
  cluster_labels[as.character(E_cluster)] <- "E"
  cluster_labels[as.character(M_cluster)] <- "M"
  if (length(other_clusters) > 0) {
    cluster_labels[as.character(other_clusters)] <- paste0("EM", seq_along(other_clusters))
  }
  
  # Map labels back to samples
  final_labels <- cluster_labels[as.character(cl1)]
  names(final_labels) <- names(cl1)
  
  return(final_labels)
}

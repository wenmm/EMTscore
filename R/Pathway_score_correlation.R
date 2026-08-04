#' Calculate correlations between pathway or signature scores across samples
#'
#' @description
#' This function calculates pairwise correlations between all numeric columns
#' of two input dataframes that contain sample-level scores (e.g., pathway, signature,
#' or module scores). The correlation is computed only for samples (rows) shared
#' by both dataframes.
#'
#' @param score_mat1 A dataframe or matrix of sample-level scores.
#' Each row represents a sample, and each column represents a pathway or signature.
#' @param score_mat2 Another dataframe or matrix of sample-level scores with matching sample row names.
#' @param method Correlation method, one of \code{"pearson"} (default), \code{"spearman"}, or \code{"kendall"}.
#'
#' @return A dataframe containing:
#' \describe{
#'   \item{Pathway_in_score_mat1}{Column name from \code{score_mat1}}
#'   \item{Pathway_in_score_mat2}{Column name from \code{score_mat2}}
#'   \item{Correlation}{Correlation coefficient}
#'   \item{P_value}{P-value of the correlation test}
#'   \item{Valid_n}{Number of non-missing paired samples used in the test}
#' }
#'
#' @details
#' This function is designed for comparing pathway, EMT, or signature scores across samples.
#' It finds the intersection of sample names (row names) between the two dataframes
#' and performs correlation analysis between every pair of numeric score columns.
#'
#' @examples
#' set.seed(1)
#' ids <- paste0("sample", seq_len(20))
#' score_mat1 <- data.frame(PathwayA = rnorm(20), PathwayB = rnorm(20),
#'   row.names = ids)
#' score_mat2 <- data.frame(EMT = rnorm(20), row.names = ids)
#' correlate_sample_scores(score_mat1, score_mat2, method = "pearson")
#' @export


correlate_sample_scores <- function(score_mat1, score_mat2, method = "pearson") {
  # Validate correlation method
  method <- match.arg(method, c("pearson", "spearman", "kendall"))
  
  # Identify common samples
  common_samples <- intersect(rownames(score_mat1), rownames(score_mat2))
  if (length(common_samples) == 0) {
    stop("No common samples found between score_mat1 and score_mat2 (row names).")
  }
  
  # Subset to common samples
  score_mat1 <- as.data.frame(score_mat1[common_samples, , drop = FALSE])
  score_mat2 <- as.data.frame(score_mat2[common_samples, , drop = FALSE])
  
  # Initialize list to store correlation results
  result_list <- list()
  
  # Compute pairwise correlations
  for (col1 in colnames(score_mat1)) {
    x <- score_mat1[[col1]]
    if (!is.numeric(x)) next
    
    for (col2 in colnames(score_mat2)) {
      y <- score_mat2[[col2]]
      if (!is.numeric(y)) next
      
      valid_idx <- complete.cases(x, y)
      if (sum(valid_idx) > 2) {
        test <- tryCatch(
          cor.test(x[valid_idx], y[valid_idx], method = method),
          error = function(e) {
            message("Invalid correlation for ", col1, " vs ", col2, ": ", e$message)
            return(NULL)
          }
        )
        
        if (!is.null(test)) {
          result_list[[length(result_list) + 1]] <- data.frame(
            Pathway_in_score_mat1 = col1,
            Pathway_in_score_mat2 = col2,
            Correlation = as.numeric(test$estimate),
            P_value = test$p.value
          )
        }
      }
    }
  }
  
  if (length(result_list) == 0) {
    stop("No valid correlations computed. Check if numeric scores exist and have overlapping samples.")
  }
  
  # Combine all results into a single dataframe
  correlation_results <- do.call(rbind, result_list)
  
  # Sort: strongest absolute correlation first, then smallest p-value
  correlation_results <- correlation_results[order(-abs(correlation_results$Correlation), correlation_results$P_value), ]
  
  rownames(correlation_results) <- NULL
  return(correlation_results)
}
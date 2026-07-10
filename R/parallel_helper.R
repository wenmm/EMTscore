#' Score every gene set of a collection in parallel
#'
#' Internal helper that applies a per-gene-set scoring function over the gene
#' sets of a collection using \pkg{BiocParallel} and row-binds the results.
#' Gene sets whose scoring function returns \code{NULL} (for example a size
#' mismatch) are dropped.
#'
#' @param n_sets Integer. Number of gene sets to iterate over.
#' @param score_one Function of a single integer index \code{k} returning a
#'   \code{data.frame} of scores for gene set \code{k}, or \code{NULL}.
#' @param cores Integer. Number of workers for
#'   \code{BiocParallel::MulticoreParam()}.
#'
#' @return A \code{data.frame} of the row-bound per-gene-set results.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @keywords internal
#' @noRd
bind_genesets_parallel <- function(n_sets, score_one, cores = 1L) {
  bpparam <- BiocParallel::MulticoreParam(workers = max(1L, as.integer(cores)))
  res <- BiocParallel::bplapply(seq_len(n_sets), score_one, BPPARAM = bpparam)
  res <- res[!vapply(res, is.null, logical(1))]
  do.call(rbind, res)
}

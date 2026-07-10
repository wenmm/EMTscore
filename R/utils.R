#' Read GMT File and Extract Genes
#'
#' This function reads a GMT (Gene Matrix Transposed) file with
#' \code{GSEABase::getGmt()} and extracts all unique genes across all gene sets.
#'
#' @param fname Character string. Path to the GMT file.
#'
#' @return A data frame with a single column \code{gene}, containing all unique
#'   genes.
#'
#' @examples
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#' genes <- read_gmt(gmt_file)
#'
#' @importFrom GSEABase getGmt geneIds
#' @export
read_gmt <- function(fname) {
  gsc <- GSEABase::getGmt(fname)
  genes <- unique(unlist(GSEABase::geneIds(gsc), use.names = FALSE))
  data.frame(gene = genes, stringsAsFactors = FALSE)
}

#' Filter a GMT file by overlap with a reference GMT
#'
#' Removes (or keeps) gene sets from `target_gmt` according to the proportion
#' of their genes that appear anywhere in `ref_gmt`.
#' By default, only novel and non-redundant gene sets are kept
#' (overlap < `cutoff`).
#'
#' @param ref_gmt Path to the reference GMT file containing the gene universe
#'   to compare against.
#' @param target_gmt Path to the GMT file that should be filtered.
#' @param output_gmt Path where the filtered GMT will be written.
#' @param cutoff Numeric threshold between 0 and 1 (default `0.5`).
#'   Fraction of overlapping genes.
#' @param keep_low_overlap Logical (default `TRUE`):
#'   \itemize{
#'     \item `TRUE`  keep gene sets with overlap < `cutoff` (remove redundant sets)
#'     \item `FALSE` keep gene sets with overlap >= `cutoff` (keep similar sets)
#'   }
#' @param min_genes Integer. Minimum number of genes required in a set
#'   (default `5`). Sets with fewer genes are ignored.
#'
#' @return Invisibly a `data.frame` with one row per processed gene set and
#'   the columns:
#'   \item{name}{Gene set name}
#'   \item{size}{Total number of genes in the set}
#'   \item{overlap}{Number of genes present in the reference GMT}
#'   \item{fraction}{`overlap / size`}
#'
#' @importFrom GSEABase getGmt geneIds toGmt setName
#' @export
#'
#' @examples
#' ref_gmt <- system.file("extdata", "TianLab_collected_EMT_signatures.gmt",
#'   package = "EMTscore")
#' target_gmt <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt",
#'   package = "EMTscore")
#' output_gmt <- tempfile(fileext = ".gmt")
#' # Keep only gene sets that share < 50 % genes with the TianLab EMT collection
#' filter_gmt_by_reference(
#'   ref_gmt,
#'   target_gmt,
#'   output_gmt,
#'   cutoff = 0.50,
#'   keep_low_overlap = TRUE
#' )
filter_gmt_by_reference <- function(ref_gmt,
                                    target_gmt,
                                    output_gmt,
                                    cutoff = 0.5,
                                    keep_low_overlap = TRUE,
                                    min_genes = 5) {
  stopifnot(
    is.numeric(cutoff), cutoff >= 0, cutoff <= 1,
    is.logical(keep_low_overlap),
    is.numeric(min_genes), min_genes >= 0
  )

  ref_gsc <- GSEABase::getGmt(ref_gmt)
  target_gsc <- GSEABase::getGmt(target_gmt)
  ref_genes <- unique(unlist(GSEABase::geneIds(ref_gsc), use.names = FALSE))
  message("Unique genes in reference GMT: ", length(ref_genes))

  target_ids <- GSEABase::geneIds(target_gsc)
  sizes <- lengths(target_ids)
  overlaps <- vapply(target_ids, function(g) sum(g %in% ref_genes), integer(1))
  fractions <- ifelse(sizes > 0, overlaps / sizes, NA_real_)

  keep <- sizes >= min_genes &
    if (keep_low_overlap) fractions < cutoff else fractions >= cutoff
  keep[is.na(keep)] <- FALSE

  stats <- data.frame(
    name = vapply(target_gsc, GSEABase::setName, character(1)),
    size = sizes,
    overlap = overlaps,
    fraction = fractions,
    stringsAsFactors = FALSE
  )
  stats <- stats[sizes >= min_genes, , drop = FALSE]

  if (nrow(stats) == 0) {
    message("No valid gene sets found in target GMT!")
    return(invisible(NULL))
  }

  GSEABase::toGmt(target_gsc[keep], output_gmt)
  message("Done: ", sum(keep), " gene sets saved to ", output_gmt)

  rownames(stats) <- NULL
  invisible(stats)
}

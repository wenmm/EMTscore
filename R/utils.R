#' Read GMT File and Extract Genes
#'
#' This function reads a GMT (Gene Matrix Transposed) file and extracts all unique genes
#' across all gene sets.
#'
#' @param fname Character string. Path to the GMT file.
#'
#' @return A data frame with a single column \code{gene}, containing all unique genes.
#'
#' @examples
#' gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
#'
#' genes <- read_gmt(gmt_file)
#'
#' @export

read_gmt <- function(fname){
  gmt_lines <- readLines(fname)
  gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, split="\t")))
  gmt_genes <- lapply(gmt_list, function(x) x[3:length(x)])
  genes <- unique(unlist(gmt_genes))
  return(data.frame(gene = genes, stringsAsFactors = FALSE))
}

#' Write a Gene Set Collection to a GMT File
#'
#' This function exports a list of gene sets into a GMT (Gene Matrix Transposed) file format.  
#' Each line of the GMT file corresponds to one gene set, with the structure:
#' By default, the gene set name is used as both the name and description fields.
#'
#' @param gene_sets A named list where each element is a character vector of gene identifiers.
#'   The list names correspond to gene set names.
#' @param file A character string specifying the path to the output GMT file.
#'
#' @details
#' The GMT format is widely used by tools such as GSEA and Enrichr for representing gene set collections.
#' Each line contains the set name, a description (often a URL or repeated name), and one or more gene identifiers separated by tabs.
#'
#' @return
#' This function writes a GMT file to disk and returns \code{NULL} invisibly.
#'
#' @examples
#' gene_sets <- list(
#'   E = c("CDH1", "VIM", "SNAI1", "ZEB1"),
#'   M = c("CDK1", "CCNB1", "CDC20")
#' )
#' write_gmt(gene_sets, file = "EM_signature.gmt")
#'
#' @export
write_gmt <- function(gene_sets, file) {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
  }
  # Construct GMT lines: "name description gene1 gene2..."
  gmt_lines <- lapply(names(gene_sets), function(name) {
    genes <- gene_sets[[name]]
    if (!is.character(genes)) {
    }
    # Use name as both set name and description (second field)
    line <- paste(c(name, name, genes), collapse = "\t")
    return(line)
  })
  # Write to file
  writeLines(unlist(gmt_lines), con = file, useBytes = TRUE)
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
#'     \item `TRUE`  – keep gene sets with overlap < `cutoff` (remove redundant sets)
#'     \item `FALSE` – keep gene sets with overlap >= `cutoff` (keep similar sets)
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
#' @export
#'
#' @examples
#' ref_gmt <- system.file("extdata", "TianLab_collected_EMT_signatures.gmt", package = "EMTscore")
#' target_gmt <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
#' output_gmt = "filtered_target.gmt"
#' # Keep only gene sets that share < 30 % genes with TianLab EMT collection
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
  
  message("=== GMT Filtering – Bulletproof Edition ===\n")
  
  # ---- Safe reading function (handles Windows/Mac/UTF-8 issues) ----
  safe_read_gmt <- function(path) {
    message("Reading → ", path)
    if (!file.exists(path)) stop("File not found: ", path)
    
    lines <- tryCatch({
      readLines(path, encoding = "UTF-8", warn = FALSE)
    }, error = function(e) {
      message("  UTF-8 failed, trying latin1...")
      readLines(path, encoding = "latin1", warn = FALSE)
    })
    message("  Success: ", length(lines), " lines loaded\n")
    return(lines)
  }
  
  ref_lines    <- safe_read_gmt(ref_gmt)
  target_lines <- safe_read_gmt(target_gmt)
  
  # ---- Build reference gene universe ----
  ref_genes <- character()
  for (line in ref_lines) {
    line <- trimws(line)
    if (line == "" || grepl("^#", line)) next
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) > 2) {
      ref_genes <- c(ref_genes, parts[-(1:2)])
    }
  }
  ref_genes <- unique(ref_genes)
  message("Unique genes in reference GMT: ", length(ref_genes), "\n")
  
  # ---- Process target GMT ----
  kept <- character()
  stats <- data.frame(name = character(), size = integer(), 
                      overlap = integer(), fraction = numeric(), 
                      stringsAsFactors = FALSE)
  
  for (line in target_lines) {
    line <- trimws(line)
    if (line == "" || grepl("^#", line)) next
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < 3) next
    
    gs_name  <- parts[1]
    gs_genes <- parts[-(1:2)]
    gs_genes <- gs_genes[gs_genes != "" & !is.na(gs_genes)]
    
    if (length(gs_genes) < min_genes) next
    
    overlap_n    <- sum(gs_genes %in% ref_genes)
    overlap_frac <- overlap_n / length(gs_genes)
    
    stats <- rbind(stats, data.frame(name = gs_name,
                                     size = length(gs_genes),
                                     overlap = overlap_n,
                                     fraction = overlap_frac))
    
    if (keep_low_overlap) {
      if (overlap_frac < cutoff) kept <- c(kept, line)
    } else {
      if (overlap_frac >= cutoff) kept <- c(kept, line)
    }
  }
  
  # ---- Summary ----
  if (nrow(stats) == 0) {
    message("No valid gene sets found in target GMT!")
    return(invisible(NULL))
  }
  
  message("Overlap fraction distribution:")
  print(summary(stats$fraction))
  
  message("\nTop 10 most similar to reference (highest overlap):")
  print(head(stats[order(stats$fraction, decreasing = TRUE), ], 10))
  
  # ---- Write output ----
  writeLines(kept, output_gmt)
  message("\nDone! ", length(kept), " gene sets saved to: ", output_gmt)
  
  if (length(kept) == 0) {
    message("\nTip: Try a lower cutoff (e.g. 0.2) or set keep_low_overlap = FALSE")
  }
  
  invisible(stats)
}
#' Example gene expression matrix
#'
#' A small log2-scale bulk gene expression matrix used to demonstrate the
#' scoring functions in runnable examples. It is a subset of the CCLE / cell
#' line panel distributed with the EMTscore study, restricted to the genes that
#' appear in the bundled EMT signature gene sets (\code{EM_signature.gmt} and
#' \code{test.gmt}) so that the examples run quickly without a network download.
#'
#' @docType data
#' @usage data(geneExp)
#' @format A numeric matrix with 544 rows (genes, row names are HGNC symbols)
#'   and 120 columns (samples/cell lines, column names match the \code{name}
#'   column of \code{\link{cell_annotation_file}}). Values are log2-transformed
#'   normalized expression.
#' @source Subset of the full expression matrix deposited at Zenodo,
#'   \doi{10.5281/zenodo.19487376} (file \code{geneExp.rda}). The processing
#'   script that produced this subset is in
#'   \code{system.file("script", "make_geneExp.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(geneExp)
#' dim(geneExp)
"geneExp"


#' Panchy et al. mesenchymal signature
#'
#' Mesenchymal (M) marker gene set curated from Panchy et al., provided as an
#' example EMT signature.
#'
#' @docType data
#' @usage data(Panchy_et_al_M_signature)
#' @format A character vector of mesenchymal-state gene symbols.
#' @source Mesenchymal marker gene list curated from Panchy et al. Gene symbols
#'   were extracted from the published EMT signature; no further processing was
#'   applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Panchy_et_al_M_signature)
"Panchy_et_al_M_signature"


#' Panchy et al. epithelial signature
#'
#' Epithelial (E) marker gene set curated from Panchy et al., provided as an
#' example EMT signature.
#'
#' @docType data
#' @usage data(Panchy_et_al_E_signature)
#' @format A character vector of epithelial-state gene symbols.
#' @source Epithelial marker gene list curated from Panchy et al. Gene symbols
#'   were extracted from the published EMT signature; no further processing was
#'   applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Panchy_et_al_E_signature)
"Panchy_et_al_E_signature"


#' Cell-line annotation table
#'
#' Sample-level annotation for the cell lines in \code{\link{geneExp}}, mapping
#' each sample name to its data source and epithelial/mesenchymal cell-type
#' label. Used to demonstrate the plotting helpers.
#'
#' @docType data
#' @usage data(cell_annotation_file)
#' @format A data.frame with 118 rows and 3 columns: \code{name} (sample
#'   identifier, matching the columns of \code{geneExp}), \code{source} (data
#'   source, e.g. CCLE), and \code{celltype_annotation} (E/M state label).
#' @source Sample-level annotation accompanying the cell-line expression matrix
#'   deposited at Zenodo, \doi{10.5281/zenodo.19487376}, subset to the samples
#'   retained in \code{\link{geneExp}}.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(cell_annotation_file)
"cell_annotation_file"


#' nnPCA_Result_EMT
#'
#' This dataset contains EMT score computed using the nnPCA method.
#'
#' @docType data
#' @usage data(nnPCA_Result_EMT)
#' @format A data frame or matrix of nnPCA-based EMT score.
#' \describe{
#'   EMT score results generated using the nnPCA method.
#' }
#' @source Precomputed example output produced by running
#'   \code{\link{Execute_nnPCA}} on the bundled \code{\link{geneExp}} matrix
#'   with the \code{EM_signature.gmt} gene set shipped in
#'   \code{system.file("extdata", package = "EMTscore")}.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(nnPCA_Result_EMT)
#'
"nnPCA_Result_EMT"


#' Gene Ontology EMT-related gene sets
#'
#' A curated collection of Gene Ontology (GO) gene sets used as input for
#' EMT scoring.
#'
#' @docType data
#' @usage data(GO)
#' @format A list of gene sets, each a character vector of gene symbols.
#' @source EMT-related gene sets derived from the Gene Ontology (GO) database
#'   (\url{http://geneontology.org}).
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(GO)
"GO"


#' MSigDB Hallmark gene sets
#'
#' Hallmark gene sets from the Molecular Signatures Database (MSigDB),
#' provided as example pathway input.
#'
#' @docType data
#' @usage data(MSigDB_Hallmark)
#' @format A list of Hallmark gene sets, each a character vector of gene symbols.
#' @source Hallmark (H) gene set collection from the Molecular Signatures
#'   Database (MSigDB, \url{https://www.gsea-msigdb.org/gsea/msigdb}).
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(MSigDB_Hallmark)
"MSigDB_Hallmark"


#' Example SCSE EMT scores
#'
#' Example EMT scores computed with the Single-Cell Signature Explorer (SCSE)
#' method.
#'
#' @docType data
#' @usage data(SCSE_Result_EMT)
#' @format A data frame of SCSE-based EMT scores (rows = samples).
#' @source Precomputed example output produced by running
#'   \code{\link{Execute_SCSE}} on the bundled \code{\link{geneExp}} matrix
#'   with the \code{EM_signature.gmt} gene set shipped in
#'   \code{system.file("extdata", package = "EMTscore")}.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(SCSE_Result_EMT)
"SCSE_Result_EMT"


#' Tan et al. cell-line epithelial signature
#'
#' Epithelial (E) signature gene set derived from cell-line data in Tan et al.
#'
#' @docType data
#' @usage data(Tan_et_al_cell_line_E_signature)
#' @format A character vector of epithelial signature gene symbols.
#' @source Cell-line epithelial EMT signature curated from Tan et al. Gene
#'   symbols were extracted from the published signature; no further processing
#'   was applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Tan_et_al_cell_line_E_signature)
"Tan_et_al_cell_line_E_signature"


#' Tan et al. cell-line mesenchymal signature
#'
#' Mesenchymal (M) signature gene set derived from cell-line data in Tan et al.
#'
#' @docType data
#' @usage data(Tan_et_al_cell_line_M_signature)
#' @format A character vector of mesenchymal signature gene symbols.
#' @source Cell-line mesenchymal EMT signature curated from Tan et al. Gene
#'   symbols were extracted from the published signature; no further processing
#'   was applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Tan_et_al_cell_line_M_signature)
"Tan_et_al_cell_line_M_signature"


#' Tan et al. tumor epithelial signature
#'
#' Epithelial (E) signature gene set derived from tumor data in Tan et al.
#'
#' @docType data
#' @usage data(Tan_et_al_tumor_E_signature)
#' @format A character vector of epithelial signature gene symbols.
#' @source Tumor epithelial EMT signature curated from Tan et al. Gene symbols
#'   were extracted from the published signature; no further processing was
#'   applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Tan_et_al_tumor_E_signature)
"Tan_et_al_tumor_E_signature"


#' Tan et al. tumor mesenchymal signature
#'
#' Mesenchymal (M) signature gene set derived from tumor data in Tan et al.
#'
#' @docType data
#' @usage data(Tan_et_al_tumor_M_signature)
#' @format A character vector of mesenchymal signature gene symbols.
#' @source Tumor mesenchymal EMT signature curated from Tan et al. Gene symbols
#'   were extracted from the published signature; no further processing was
#'   applied.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(Tan_et_al_tumor_M_signature)
"Tan_et_al_tumor_M_signature"


#' Example data prepared for plotting
#'
#' A data frame combining cell annotations with computed scores, ready to be
#' passed to the plotting helpers.
#'
#' @docType data
#' @usage data(data_for_plot)
#' @format A data frame (rows = cells/samples) of annotations and scores.
#' @source Precomputed example object built with \code{\link{data_prepare}} by
#'   merging \code{\link{cell_annotation_file}} with nnPCA scores computed from
#'   the bundled \code{\link{geneExp}} matrix.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(data_for_plot)
"data_for_plot"


#' Example multiple nnPCA results
#'
#' Example Escore and Mscore values computed with the nnPCA method across
#' multiple gene sets.
#'
#' @docType data
#' @usage data(nnPCA_Result_multiple)
#' @format A data frame with nnPCA Escore and Mscore columns (rows = samples).
#' @source Precomputed example output produced by running \code{\link{Execute_nnPCA}}
#'   (with \code{dimension = 2}) on the bundled \code{\link{geneExp}} matrix and
#'   the \code{EM_signature.gmt} gene set shipped in
#'   \code{system.file("extdata", package = "EMTscore")}.
#'   Creation and provenance code: \code{system.file("script", "make-data.R", package = "EMTscore")}.
#' @keywords datasets
#' @examples
#' data(nnPCA_Result_multiple)
"nnPCA_Result_multiple"

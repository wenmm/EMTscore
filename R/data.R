

#' Panchy_et_al_M_signature
#'
#' This is an example dataset.
#'
#' @docType data
#' @usage data(Panchy_et_al_M_signature)
#' @format M gene list from Panchy_et_al
#' \describe{
#'   M gene
#' }
#' @examples
#' data(Panchy_et_al_M_signature)
#' 
"Panchy_et_al_M_signature"


#' Panchy_et_al_E_signature
#'
#' This is an example dataset.
#'
#' @docType data
#' @usage data(Panchy_et_al_E_signature)
#' @format E gene list from Panchy_et_al
#' \describe{
#'   E gene
#' }
#' @examples
#' data(Panchy_et_al_E_signature)
#' 
"Panchy_et_al_E_signature"


#' cell_annotation_file
#'
#' This is a cell type annotation file
#'
#' @docType data
#' @usage data(cell_annotation_file)
#' @format example data
#' \describe{
#'  celltype annotation
#' }
#' @examples
#' data(cell_annotation_file)
#' 
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
#' @keywords datasets
#' @examples
#' data(nnPCA_Result_multiple)
"nnPCA_Result_multiple"

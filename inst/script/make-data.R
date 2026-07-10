## Provenance and creation code for the datasets in EMTscore's data/ directory.
##
## This script documents where every object under data/ comes from and, where
## the object is derived from other bundled data, gives runnable code that
## reproduces it exactly. It is not run at build/install time; it is included
## for provenance per the Bioconductor package guidelines.
##
## Objects fall into three groups:
##   (A) Curated gene signatures / gene sets obtained from external sources.
##   (B) Sample data (expression matrix and its annotation) from the study's
##       Zenodo deposit.
##   (C) Example score objects that are DERIVED by running EMTscore's own
##       functions on the bundled data (fully reproducible below).

library(EMTscore)


## ---------------------------------------------------------------------------
## (A) Curated gene signatures / gene sets (external sources)
## ---------------------------------------------------------------------------
## Each is stored as a one-column data.frame with column `GeneName`. The gene
## symbols were taken from the published signatures / databases cited below; no
## further processing beyond extracting the symbol list was applied.
##
##   Panchy_et_al_E_signature, Panchy_et_al_M_signature
##       Epithelial / mesenchymal marker gene lists from Panchy et al.
##       (see the EMTscore vignette reference "Panchy et al.").
##       License: as per the source publication.
##
##   Tan_et_al_cell_line_E_signature, Tan_et_al_cell_line_M_signature,
##   Tan_et_al_tumor_E_signature,     Tan_et_al_tumor_M_signature
##       Generic epithelial / mesenchymal EMT signatures (cell-line and tumor
##       versions) from Tan et al. 2014, EMBO Molecular Medicine.
##       License: as per the source publication.
##
##   MSigDB_Hallmark
##       The HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION gene set (200 genes)
##       from the MSigDB Hallmark collection (Liberzon et al. 2015). It is the
##       gene list of the Hallmark EMT .gmt shipped in inst/extdata and is
##       reproduced exactly from that file:
##   MSigDB_Hallmark <- data.frame(
##     GeneName = read_gmt(system.file(
##       "extdata", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.gmt",
##       package = "EMTscore"))$gene,
##     stringsAsFactors = FALSE)
##
##   GO
##       Epithelial-mesenchymal-transition-related genes from the Gene Ontology
##       (Gene Ontology Consortium; http://geneontology.org), 90 gene symbols.
##       License: Gene Ontology data are distributed under CC BY 4.0.


## ---------------------------------------------------------------------------
## (B) Sample expression data and its annotation
## ---------------------------------------------------------------------------
##   geneExp
##       A 544 x 120 subset of the full CCLE / cell-line expression matrix
##       deposited at Zenodo (doi:10.5281/zenodo.19487376, file geneExp.rda),
##       restricted to the genes of the bundled EMT signature .gmt files so the
##       examples run quickly. See inst/script/make_geneExp.R for the exact
##       subsetting code.
##
##   cell_annotation_file
##       Sample-level annotation (columns name / source / celltype_annotation)
##       accompanying the expression matrix at the same Zenodo deposit
##       (doi:10.5281/zenodo.19487376), subset to the samples kept in geneExp.


## ---------------------------------------------------------------------------
## (C) Derived example score objects (reproducible from the bundled data)
## ---------------------------------------------------------------------------
## These are precomputed outputs of EMTscore's scoring/plotting functions on
## the bundled geneExp matrix, stored so that the plotting examples run without
## recomputation. The code below regenerates each object exactly (verified to
## match the shipped objects).

data(geneExp, package = "EMTscore")
data(cell_annotation_file, package = "EMTscore")

hallmark_gmt <- system.file(
  "extdata", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.gmt",
  package = "EMTscore")
em_gmt <- system.file("extdata", "EM_signature.gmt", package = "EMTscore")

## Single-gene-set EMT scores (one column "EMT", rows = samples)
nnPCA_Result_EMT <- Execute_nnPCA(geneExp, hallmark_gmt,
  dimension = 1, score_names = "EMT")
SCSE_Result_EMT <- Execute_SCSE(geneExp, hallmark_gmt, score_names = "EMT")

## Multi-gene-set nnPCA scores (samples x gene sets)
nnPCA_Result_multiple <- Execute_nnPCA_parallel(geneExp, em_gmt,
  dimension = 1, cores = 1)

## Annotation merged with the multi-gene-set scores, ready for the plot helpers
data_for_plot <- data_prepare(cell_annotation_file, nnPCA_Result_multiple,
  merge_colname = "name")

## Objects are saved with usethis::use_data(..., compress = "xz") /
## save(..., compress = "xz") into data/.

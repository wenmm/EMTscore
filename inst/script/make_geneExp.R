# Provenance script for data/geneExp.rda
#
# The example expression matrix `geneExp` is a subset of the full CCLE / cell
# line expression matrix deposited at Zenodo (doi:10.5281/zenodo.19487376,
# file `geneExp.rda`, 15950 genes x 120 samples, log2-normalized expression).
#
# To keep the shipped example small and to let the man-page examples run
# without a network download, the matrix is restricted to the genes contained
# in the two bundled EMT signature gene sets (EM_signature.gmt and test.gmt).
#
# Reproduce with:

url <- "https://zenodo.org/records/19487376/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile) # loads `geneExp` (15950 x 120)
full <- geneExp

gmt_genes <- function(f) {
  path <- system.file("extdata", f, package = "EMTscore")
  lines <- readLines(path)
  unique(unlist(lapply(lines, function(x) {
    fields <- strsplit(x, "\t")[[1]]
    fields[-(1:2)]
  })))
}

genes <- unique(c(gmt_genes("EM_signature.gmt"), gmt_genes("test.gmt")))
genes <- intersect(genes, rownames(full))

geneExp <- full[genes, , drop = FALSE] # 544 x 120
save(geneExp, file = "data/geneExp.rda", compress = "xz")

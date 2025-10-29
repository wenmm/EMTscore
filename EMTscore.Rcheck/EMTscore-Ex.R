pkgname <- "EMTscore"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "EMTscore-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('EMTscore')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Execute_AUCell")
### * Execute_AUCell

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_AUCell
### Title: Execute AUCell Scoring from GMT File
### Aliases: Execute_AUCell

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
result <- Execute_AUCell(geneExp, gmt_file, score_names = "score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_AUCell", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_AUCell_parallel")
### * Execute_AUCell_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_AUCell_parallel
### Title: Compute AUCell Scores for Multiple Gene Sets
### Aliases: Execute_AUCell_parallel

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
result <- Execute_AUCell_parallel(geneExp, gmt_file, cores = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_AUCell_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_JAS")
### * Execute_JAS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_JAS
### Title: Execute JASMINE Scoring from GMT File
### Aliases: Execute_JAS

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
result <- Execute_JAS(geneExp, gmt_file, score_names = "score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_JAS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_JASMINE_parallel")
### * Execute_JASMINE_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_JASMINE_parallel
### Title: Compute JASMINE Scores for Multiple Gene Sets in Parallel
### Aliases: Execute_JASMINE_parallel

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
result <- Execute_JASMINE_parallel(geneExp, gmt_file, cores = 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_JASMINE_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_SCSE")
### * Execute_SCSE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_SCSE
### Title: Execute SCSE Scoring from GMT File
### Aliases: Execute_SCSE

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
result <- Execute_SCSE(geneExp, gmt_file, score_names = "score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_SCSE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_SCSE_parallel")
### * Execute_SCSE_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_SCSE_parallel
### Title: Compute SCSE Scores for Multiple Gene Sets
### Aliases: Execute_SCSE_parallel

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
result <- Execute_SCSE_parallel(data, gmt_file, cores = 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_SCSE_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_nnPCA")
### * Execute_nnPCA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_nnPCA
### Title: Execute Non-Negative Sparse PCA (nnPCA) Scoring
### Aliases: Execute_nnPCA

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
result <- Execute_nnPCA(geneExp, gmt_file,dimension = 1, score_names = "score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_nnPCA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_nnPCA_parallel")
### * Execute_nnPCA_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_nnPCA_parallel
### Title: Compute nnPCA Scores for Multiple Gene Sets in Parallel
### Aliases: Execute_nnPCA_parallel

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
result <- Execute_nnPCA_parallel(geneExp, gmt_file, dimension = 1, cores = 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_nnPCA_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_ssGSEA_parallel")
### * Execute_ssGSEA_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_ssGSEA_parallel
### Title: Compute ssGSEA Scores for Multiple Gene Sets in Parallel
### Aliases: Execute_ssGSEA_parallel

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "h.all.v2025.1.Hs.symbols.gmt", package = "EMTscore")
result <- Execute_ssGSEA_parallel(geneExp, gmt_file, cores = 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_ssGSEA_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Execute_ssGSVA")
### * Execute_ssGSVA

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Execute_ssGSVA
### Title: Execute ssGSVA Scoring from GMT File
### Aliases: Execute_ssGSVA

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
result <- Execute_ssGSVA(geneExp, gmt_file, score_names = "score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Execute_ssGSVA", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("add_EMT_score")
### * add_EMT_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: add_EMT_score
### Title: Add EMT Scores to Seurat Objects
### Aliases: add_EMT_score

### ** Examples

rds_file <- system.file("extdata", "example_seurat_obj.rds", package = "EMTscore")
seurat_files <- c("/home/hwen6/project/EMTscore/2025_09_24/EMTscore/inst/extdata/example_seurat_obj.rds", "/home/hwen6/project/EMTscore/2025_09_24/EMTscore/inst/extdata/example_seurat_obj.rds")
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

# Using nnPCA
seurat_list <- add_EMT_score(seurat_files, gmt_file, emt_name = "EMT_score", method = "nnPCA")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("add_EMT_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("add_EMT_score_Bulk")
### * add_EMT_score_Bulk

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: add_EMT_score_Bulk
### Title: Compute EMT scores for one or multiple bulk expression matrices
### Aliases: add_EMT_score_Bulk

### ** Examples

# single matrix
url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")
res <- add_EMT_score_Bulk(expr_mat = geneExp, gmt_file, method = "AUCell")

# multiple matrices (named list)
exprs <- list(TCGA = geneExp, GTEX = geneExp)
res_list <- add_EMT_score_Bulk(expr_mat_list = exprs, gmt_file, method = "nnPCA", dimension = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("add_EMT_score_Bulk", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_Signature_score")
### * compute_Signature_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_Signature_score
### Title: Compute Signature Score from Bulk Matrix or Seurat Object
### Aliases: compute_Signature_score

### ** Examples

url <- "https://zenodo.org/records/17438655/files/geneExp.rda"
destfile <- tempfile(fileext = ".rda")
download.file(url, destfile, mode = "wb")
load(destfile)
signature_file <- system.file("extdata", "stemsig.tsv", package = "EMTscore")
compute_Signature_score(geneExp, signature_file, "Stemness_Score")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_Signature_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("correlate_sample_scores")
### * correlate_sample_scores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: correlate_sample_scores
### Title: Calculate correlations between pathway or signature scores
###   across samples
### Aliases: correlate_sample_scores

### ** Examples

# score_mat1: pathway scores per sample; score_mat2: EMT signature scores per sample
data(SCSE_Result_EMT)
data(SCSE_Result_multiple)
result <- correlate_sample_scores(score_mat1 = SCSE_Result_multiple, score_mat2 = SCSE_Result_EMT, method = "spearman")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("correlate_sample_scores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_EMT_from_objects")
### * plot_EMT_from_objects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_EMT_from_objects
### Title: Plot EMT Scores Along Pseudotime from Seurat Objects
### Aliases: plot_EMT_from_objects

### ** Examples

# Assume seurat_list has EMT scores calculated
p <- plot_EMT_from_objects(seurat_list, col_name = "Pseudotime", emt_score_col = "EMT_Score")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_EMT_from_objects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read_gmt")
### * read_gmt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read_gmt
### Title: Read GMT File and Extract Genes
### Aliases: read_gmt

### ** Examples

gmt_file <- system.file("extdata", "test.gmt", package = "EMTscore")

genes <- read_gmt(gmt_file)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read_gmt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("write_gmt")
### * write_gmt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: write_gmt
### Title: Write Gene Sets to a GMT File
### Aliases: write_gmt

### ** Examples

gene_sets <- list(
  EMT_High = c("VIM", "ZEB1", "SNAI1"),
  EMT_Low  = c("CDH1", "OCLN", "DSP")
)
write_gmt(gene_sets, "EMT_signature.gmt")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("write_gmt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

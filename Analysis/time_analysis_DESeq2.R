# Time course analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)

## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start 
# Use test="LRT" for significance testing when working with single-cell data, over the Wald test. 
#Set the following DESeq arguments to these values: useT=TRUE, minmu=1e-6, and minReplicatesForReplace=Inf. The default setting of minmu was benchmarked on bulk RNA-seq and is not appropriate for single cell data when the expected count is often much less than 1.
#The default size factors are not optimal for single cell count matrices, instead consider setting sizeFactors from scran::computeSumFactors.

# another useful link: https://bioconductor.org/books/3.12/OSCA/single-nuclei-rna-seq-processing.html

library(readr)
common_genes_with_avg_log2FC <- read_csv("common_genes_with_avg_log2FC.csv")
View(common_genes_with_avg_log2FC)

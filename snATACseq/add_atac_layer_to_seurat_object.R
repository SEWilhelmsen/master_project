# Add ATAC layer to seurat data
# Pre-processing of ATAC-sequencing data. 

# Load libraries
#############################################################################
library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
#library(EnsDb.Hsapiens.v86) # For human annotations
library(EnsDb.Mmusculus.v79) # For mouse annotations
library(Seurat)


## Set parameters                                                         
#################################################################################

## .h5 file to read
h5_file <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Multiome seq data/TP2-4/outs/filtered_feature_bc_matrix.h5" 

# set a project name to name saved output files and the Seurat object
sample_name <- "TP2-4" 

## Remember to create an output directory for saving files and plots
#output.dir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Multiome seq data/QC/TP2-4/"
output_dir <- "C:/Users/siljeew/snATACseq/Analysis/"

# Preparation
################################################################################

# Load data and create Seurat object
dataset <- Read10X_h5(h5_file, use.names = TRUE, unique.features = TRUE)

# extract RNA and ATAC data
rna_counts <- dataset$`Gene Expression`
atac_counts <- dataset$Peaks

# Create seurat object to store RNA counts and metadata
SeuratObj <- CreateSeuratObject(counts = rna_counts, project = sample_name)

# Add data to the seurat object
###############################################################################
# add percent mitochondrial genes
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^mt-")
# add percent ribosomal genes
SeuratObj[["percent.ribo"]] <- PercentageFeatureSet(SeuratObj, pattern = "^Rp[sl]") 

# Add the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # Genomic annotation and fragment files
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# Add assay to the seurat object as another layer of data
#frag.file <- "/Volumes/Backup_AB/GCF0742_GCF0795_CellRanger/Sample1-TP75B/outs/atac_fragments.tsv.gz"
frag.file <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Multiome seq data/TP2-4/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
SeuratObj[["ATAC"]] <- chrom_assay

# Create plots
#####################################################

# Visualize QC metrics
VlnPlot(SeuratObj, features = c("nCount_ATAC", "nFeature_ATAC", "TSS.enrichment"), ncol = 3)

# Create UMAP plot
source("C:/Users/siljeew/snATACseq/Visualization/create_umap_plot.R")


# # Check assay slots
# GetAssayData(SeuratObj[['ATAC']], slot = "data")
# GetAssayData(SeuratObj[['ATAC']], slot = "counts")


# Clustering
SeuratObj <- FindNeighbors(SeuratObj, reduction = "lsi", dims = 1:30)
SeuratObj <- FindClusters(SeuratObj, resolution = 0.5)

# Feature Plot
FeaturePlot(SeuratObj, features = "A specific peak or gene")


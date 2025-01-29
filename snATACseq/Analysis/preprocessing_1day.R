### ATAC sequencing data SHAM vs. AB - Preprocessing
## Timepoint: 1 Day

library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(Seurat)


## TP17-2

h5_file <- "W:/CardioTarget/Multiome seq data/TP17-2/outs/filtered_feature_bc_matrix.h5"
sample_name <- "TP17-2"
condition_and_time <- "1 Day - SHAM"

output.dir <- "C:/Users/siljeew/snATACseq/"

# Load data and create seurat object
dataset <- Read10X_h5(h5_file, use.names = TRUE, unique.features = TRUE)

# Extract rna and atac data
rna_counts <- dataset$`Gene Expression`
atac_counts <- dataset$Peaks

# Put rna and atac data into a seurat object
dataset_seurat <- CreateSeuratObject(counts = rna_counts, project = sample_name) # Should project be time point and condition instead? 

# Add percentage of mitochondrial genes
dataset_seurat[["percent.mt"]] <- PercentageFeatureSet(dataset_seurat, pattern = "^mt-")
# Add percentage of ribosomal genes
dataset_seurat[["percent.ribo"]] <- PercentageFeatureSet(dataset_seurat, pattern = "^Rp[sl]") 

# Add the ATAC-seq data
grange_counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
atac_counts <- atac_counts[as.vector(grange_use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) ## human = EnsDb.Hsapiens.v86
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10" # Genome version mm10

print(seqlevelsStyle(annotations))
print(genome(annotations))


frag.file <- "W:/CardioTarget/Multiome seq data/TP17-2/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
dataset_seurat[["ATAC"]] <- chrom_assay

# View metadata columns
head(dataset_seurat@meta.data)

# Check for specific columns
colnames(dataset_seurat@meta.data)

# Add condition and time point
dataset_seurat$condition_and_time <- condition_and_time

# View updated metadata
head(dataset_seurat@meta.data)

# Optionally view specific columns
colnames(dataset_seurat@meta.data)

# Quality control
##########################################################
# Set correct assay
if ("ATAC" %in% Assays(dataset_seurat)) {
  DefaultAssay(dataset_seurat) <- "ATAC"
} else {
  stop("ATAC assay not found in the dataset.")
}

# Calculate TSS enrichment 
if (is(dataset_seurat[["ATAC"]], "ChromatinAssay")) {
  
  # Calculate TSS enrichment
  dataset_seurat <- TSSEnrichment(object = dataset_seurat, fast = FALSE)
} else {
  stop("The ATAC assay is not a ChromatinAssay")
}

head(dataset_seurat@meta.data) # Check if TSS is available

# Violin plot for TSS Enrichment
VlnPlot(dataset_seurat, features = "TSS.enrichment", pt.size = 0.1)

# Summary of TSS Enrichment
summary(dataset_seurat$TSS.enrichment)

# Nucleosome signal
dataset_seurat <- NucleosomeSignal(dataset_seurat)

View(dataset_seurat@meta.data)

# Investigate the data based on violin pltots
VlnPlot(object = dataset_seurat, 
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
        pt.size = 0.1,
        ncol = 4)


dataset_seurat_subset <- subset(x = dataset_seurat,
               subset = nCount_ATAC > 3000 &
                 nCount_ATAC < 30000 &
                 #pct_reads_in_peaks > 15 & 
                 #blacklist_ratio < 0.05 &
                 nucleosome_signal < 4 &
                 TSS.enrichment > 3)



# Normalization and linear dimensional reduction 
dataset_seurat_subset <- RunTFIDF(dataset_seurat_subset) # normalization
dataset_seurat_subset <- FindTopFeatures(dataset_seurat_subset, min.cutoff = 'q0') # selecting top features
dataset_seurat_subset <- RunSVD(dataset_seurat_subset) # dimensionality reduction

DepthCor(dataset_seurat_subset) # The first column represents technical something. Only use columns 2 +

# Non-linear dimensional reduction and clustering 
dataset_seurat_subset <- RunUMAP(object = dataset_seurat_subset, reduction = 'lsi', dims = 2:30)
dataset_seurat_subset <- FindNeighbors(object = dataset_seurat_subset, reduction = 'lsi', dims = 2:30)
dataset_seurat_subset <- FindClusters(object = dataset_seurat_subset, algorithm = 3)

DimPlot(object = dataset_seurat_subset, label = TRUE) + NoLegend()
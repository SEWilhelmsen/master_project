
library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(Seurat)

# Define file locations for each sample
samples <- list(
  TP172 = list (h5_file = "W:/CardioTarget/Multiome seq data/TP17-2/outs/filtered_feature_bc_matrix.h5",
                frag_file = "W:/CardioTarget/Multiome seq data/TP17-2/outs/atac_fragments.tsv.gz",
                project = "1 Day - SHAM"),
  TP173 = list(h5_file = "W:/CardioTarget/Multiome seq data/TP17-3/outs/filtered_feature_bc_matrix.h5",
               frag_file = "W:/CardioTarget/Multiome seq data/TP17-3/outs/atac_fragments.tsv.gz",
               project= "1 Day - SHAM"),
  TP131 = list(h5_file = "W:/CardioTarget/Multiome seq data/TP13-1/outs/filtered_feature_bc_matrix.h5",
               frag_file = "W:/CardioTarget/Multiome seq data/TP13-1/outs/atac_fragments.tsv.gz",
               project= "1 Day - SHAM"),
  TP142 = list(h5_file = "W:/CardioTarget/Multiome seq data/TP14-2/outs/filtered_feature_bc_matrix.h5",
               frag_file = "W:/CardioTarget/Multiome seq data/TP14-2/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB"),
  TP141 = list(h5_file = "W:/CardioTarget/Multiome seq data/TP14-1/outs/filtered_feature_bc_matrix.h5",
               frag_file = "W:/CardioTarget/Multiome seq data/TP14-1/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB"),
  TP153 = list(h5_file = "W:/CardioTarget/Multiome seq data/TP15-3/outs/filtered_feature_bc_matrix.h5",
               frag_file = "W:/CardioTarget/Multiome seq data/TP15-3/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB")
)


# Function to create a seurat object with RNA and ATAC data
seurat_list <- mapply(function(sample, sample_name) {
  # Control which sample is being processed
  message("Processing sample: ", sample_name)
  
  #Load data
  dataset <- Read10X_h5(sample$h5_file, use.names = TRUE, unique.features = TRUE)
  
  # extract RNA and ATAC data
  rna_counts <- dataset$`Gene Expression`
  atac_counts <- dataset$Peaks
  
  # Create seurat object with RNA assay
  seurat_obj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = sample$project
  )
  
  # Add ATAC as a new assay
  frag_file <- sample$frag_file
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = frag_file,
    annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
    genome = 'mm10'
  )
  
  # Add chromatin assay
  seurat_obj [["ATAC"]] <- chrom_assay 

  
  return(seurat_obj)
}, sample = samples, sample_name = names(samples), SIMPLIFY = FALSE)

# Merge all seurat objects into one
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)

print(Assays(combined_seurat))
head(combined_seurat@meta.data)

# Add data to the seurat object
# Add percent mitochondrial genes
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")

# Add percent ribosomal genes
combined_seurat[["percent.ribo"]] <- PercentageFeatureSet(combined_seurat, pattern = "^Rp[sl]") 
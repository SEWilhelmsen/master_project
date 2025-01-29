# Preprocessing ATAC data for six samples
# Time point:  1 day

library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(Seurat)



# Define file locations for each sample
samples <- list(
  TP172 = list (h5_file = "D:/MultiomeSeq/TP17-2/outs/filtered_feature_bc_matrix.h5",
                frag_file = "D:/MultiomeSeq/TP17-2/outs/atac_fragments.tsv.gz",
                project = "1 Day - SHAM"),
  TP173 = list(h5_file = "D:/MultiomeSeq/TP17-3/outs/filtered_feature_bc_matrix.h5",
               frag_file = "D:/MultiomeSeq/TP17-3/outs/atac_fragments.tsv.gz",
               project= "1 Day - SHAM"),
  TP131 = list(h5_file = "D:/MultiomeSeq/TP13-1/outs/filtered_feature_bc_matrix.h5",
               frag_file = "D:/MultiomeSeq/TP13-1/outs/atac_fragments.tsv.gz",
               project= "1 Day - SHAM"),
  TP142 = list(h5_file = "D:/MultiomeSeq/TP14-2/outs/filtered_feature_bc_matrix.h5",
               frag_file = "D:/MultiomeSeq/TP14-2/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB"),
  TP141 = list(h5_file = "D:/MultiomeSeq/TP14-1/outs/filtered_feature_bc_matrix.h5",
               frag_file = "D:/MultiomeSeq/TP14-1/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB"),
  TP153 = list(h5_file = "D:/MultiomeSeq/TP15-3/outs/filtered_feature_bc_matrix.h5",
               frag_file = "D:/MultiomeSeq/TP15-3/outs/atac_fragments.tsv.gz",
               project= "1 Day - AB")
)

# Retrieve annotations with the correct genome
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"         



# Create a seurat object with ATAC and RNA data 
##################################################################################
# Function to create a seurat object with RNA and ATAC data
seurat_list <- mapply(function(sample, sample_name) {
  # Control which sample is being processed
  message("Processing sample: ", sample_name)
  
  #Load data
  dataset <- Read10X_h5(sample$h5_file, use.names = TRUE, unique.features = TRUE)
  
  # extract RNA and ATAC data
  rna_counts <- dataset$`Gene Expression`
  atac_counts <- dataset$Peaks
  
  # Reduce noise and potential artifacts
  # Convert string names to genomic ranges
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  # Filter for standard chromosome names
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  # Filter out peaks that don't reside on standard chromosomes
  atac_counts <- atac_counts[as.vector(grange.use), ]
  
  
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
    fragments = sample$frag_file,
    annotation = annotations,
    genome = 'mm10'
  )
  
  # Add chromatin assay
  seurat_obj [["ATAC"]] <- chrom_assay 

  
  return(seurat_obj)
}, sample = samples, sample_name = names(samples), SIMPLIFY = FALSE)

# Merge all seurat objects into one. 
# Seurat appends a sample/specific prefix by default to handle duplicated names across objects, you dont have to worry.
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)

# Inspect data
print(Assays(combined_seurat))
head(combined_seurat@meta.data)



# Add meta data to the seurat object
####################################################################################
# Add percent mitochondrial genes
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^mt-")
head(combined_seurat@meta.data)

# Add percent ribosomal genes
combined_seurat[["percent.ribo"]] <- PercentageFeatureSet(combined_seurat, pattern = "^Rp[sl]") 
head(combined_seurat@meta.data)

# Add nucleosome signal score per cell - maybe do this with the separate samples and not the combined? 
combined_seurat <- NucleosomeSignal(combined_seurat, assay = "ATAC", verbose = TRUE)
head(combined_seurat@meta.data)

# Add TSS enrichment score per cell
combined_seurat <- TSSEnrichment(object = combined_seurat, assay = "ATAC", fast = FALSE)
head(combined_seurat@meta.data)



# Visualize for quality control
############################################################
colnames(combined_seurat@meta.data)

# Create density and scatter plot
a1 <- DensityScatter(combined_seurat, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(combined_seurat, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
combined_plot <- a1 | a2 + plot_annotation(title = "ATAC-seq Quality Control Metrics at 1 Day") # Print the two plots
print(combined_plot)

# Create violin plots
VlnPlot(object = combined_seurat,
        features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
        pt.size = 0.1,
        ncol = 4,
        assay = "ATAC")

# Save seurat object before filtering
SaveSeuratRds(combined_seurat, file = "D:/ATAC/combined_seurat_1d.Rds")



# Filter based on quality 
############################################################
combined_seurat_filtered <- subset(x = combined_seurat,
                                   subset = nCount_ATAC > 3000 &
                                     nCount_ATAC < 30000 &
                                     nucleosome_signal < 4 &
                                     TSS.enrichment > 3)


# List assays and set default assay
Assays(combined_seurat_filtered)
DefaultAssay(combined_seurat_filtered) <- "ATAC"

# Normalization 
combined_seurat_filtered <- RunTFIDF(combined_seurat_filtered) 

# Find top features
combined_seurat_filtered <- FindTopFeatures(combined_seurat_filtered, min.cutoff = 'q0')

# Dimensionality reduction: Singular-value decomposition
combined_seurat_filtered <- RunSVD(combined_seurat_filtered)

# Visualization
DepthCor(combined_seurat_filtered) # The first point is often correlated with sequencing depth? wtf?

# Dimensionality reduction: 
combined_seurat_filtered <- RunUMAP(object = combined_seurat_filtered, reduction = 'lsi', dims = 2:30) #Latent Semantic Indexing (LSI)
combined_seurat_filtered <- FindNeighbors(object = combined_seurat_filtered, reduction = 'lsi', dims = 2:30)

# Clustering
combined_seurat_filtered <- FindClusters(object = combined_seurat_filtered, algorithm = 3)

DimPlot(object = combined_seurat_filtered, label = TRUE) + NoLegend()


# Save seurat 
SaveSeuratRds(combined_seurat_filtered, file = "D:/ATAC/combined_seurat_filtered_1d.Rds")


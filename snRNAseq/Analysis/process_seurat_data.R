# Process seurat data 
# Silje Wilhelmsen

# Load libraries
library(dplyr)
library(tidyr)
library(Seurat)
library(ggpubr)
library(ggplot2)
library(tibble)


# Function to adjust Seurat metadata
process_metadata <- function(seurat_object) {
  seurat_object@meta.data <- seurat_object@meta.data %>%
    mutate(time_point = sub(" - .*", "", orig.ident), 
           condition = sub(".*- ", "", orig.ident),
           condition = ifelse(condition == "AB", "ORAB", condition))
  return(seurat_object)
}


# Function to extract and filter gene expression data
process_expression_data <- function(seurat_object, genes_of_interest) {
  data <- GetAssayData(seurat_object, layer = "data") %>%
    as.data.frame()
  
  # Keep only the genes from genes_of_interest
  data <- data[rownames(data) %in% toupper(genes_of_interest), ]
  print(paste("Number of genes in filtered data:", nrow(data)))
  
  rownames(data) <- toupper(rownames(data))
  
  # Convert data to long format
  data_long <- data %>%
    rownames_to_column(var = "gene") %>%
    pivot_longer(cols = -gene, names_to = "cell", values_to = "expression") %>%
    mutate(time_point = seurat_object@meta.data$time_point[match(cell, rownames(seurat_object@meta.data))],
           condition = seurat_object@meta.data$condition[match(cell, rownames(seurat_object@meta.data))])
  
  return(data_long)
}

# Process each Seurat object and return combined data
process_seurat_object <- function(file_path, genes_of_interest) {
  seurat_object <- readRDS(file_path)
  seurat_object <- process_metadata(seurat_object)
  data_long <- process_expression_data(seurat_object, genes_of_interest)
  return(data_long)
}

# Combine data from all Seurat objects
#combined_data_long <- do.call(rbind, lapply(file_paths, process_seurat_object, genes_of_interest = genes_of_interest))

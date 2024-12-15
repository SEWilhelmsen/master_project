### RNA sequencing data SHAM vs. AB - Preprocessing
## Silje Wilhelmsen

## Load the data for timepoint and store as seurat objects
create_mouse_seurat_object <- function(experiment_data, time_point) {
  seurat_object_per_mouse <- list()
  
  # Loop through each item in the list
  mouse_ids <- names(experiment_data)
  for (mouse_id in mouse_ids) {
    mouse <- experiment_data[[mouse_id]]
    path <- mouse$path
    condition <- mouse$condition
    print(paste("Mouse ID:", mouse_id, "Path:", path, "Condition:", condition))
    
    counts <- Read10X_h5(path)
    
    project <- paste(time_point, " - ", condition)
    
    mouse_seurat_object <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA", 
      project = project
     )
    
    seurat_object_per_mouse[[mouse_id]] <- mouse_seurat_object
    
  }
  
  return(seurat_object_per_mouse) 
}


## Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from
annotation_of_nuclei_mouse <- function(seurat_object_per_mouse) {
  predictions_list <- list() # Create an empty list to store the prediction results per mouse
  
  mouse_ids <- names(seurat_object_per_mouse)
  for (mouse_id in mouse_ids) {
    predictions_list[[mouse_id]] <- RunAzimuth(mouse_id, reference = "heartref")
  }
  
  return(predictions_list)
}


# Change names of rows of the predictions seurat objects to nuclei ID and animal ID
change_names_of_rows_and_columns <- function(annotation_of_nuclei_mouse(seurat_object_per_mouse)) {
  rownames_list <- list() # Create an empty list to store rownames
  
  mouse_ids <- names(seurat_object_per_mouse)
  for (mouse_id in mouse_ids) {
    
  }
  
  return(rownames_list)
  
  rownames_TP51 <- rownames(Prediction_TP51@meta.data)
  new_rownames_TP51 <- gsub("-1$", "-TP51", rownames_TP51)
  colnames(Prediction_TP51)<-new_rownames_TP51
  
}
#Prediction_TP51 <- RunAzimuth(TP51, reference = "heartref")

####################################################################

## Change names on the rows of the Predictions seurat objects to nuclei ID + Animal ID#
rownames_TP51 <- rownames(Prediction_TP51@meta.data)
rownames_TP75 <- rownames(Prediction_TP75@meta.data)
rownames_TP52 <- rownames(Prediction_TP52@meta.data)
rownames_TP511 <- rownames(Prediction_TP511@meta.data)
rownames_TP62 <- rownames(Prediction_TP62@meta.data)
rownames_TP513 <- rownames(Prediction_TP513@meta.data)

new_rownames_TP51 <- gsub("-1$", "-TP51", rownames_TP51)
new_rownames_TP75 <- gsub("-1$", "-TP75", rownames_TP75)
new_rownames_TP52 <- gsub("-1$", "-TP52", rownames_TP52)
new_rownames_TP511 <- gsub("-1$", "-TP511", rownames_TP511)
new_rownames_TP62 <- gsub("-1$", "-TP62", rownames_TP62)
new_rownames_TP513 <- gsub("-1$", "-TP513", rownames_TP513)

colnames(Prediction_TP51)<-new_rownames_TP51
colnames(Prediction_TP75)<-new_rownames_TP75
colnames(Prediction_TP52)<-new_rownames_TP52
colnames(Prediction_TP511)<-new_rownames_TP511
colnames(Prediction_TP62)<-new_rownames_TP62
colnames(Prediction_TP513)<-new_rownames_TP513

## Merge all 6 animals into one Seurat object

mouse_3w <- merge(Prediction_TP51,
                  y=c(Prediction_TP75, Prediction_TP52, Prediction_TP511, Prediction_TP62, Prediction_TP513),
                  project = "Full Merge Project")
mouse_3w <- NormalizeData(mouse_3w,
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
mouse_3w <- FindVariableFeatures(mouse_3w)
mouse_3w <- ScaleData(mouse_3w)
mouse_3w <- RunPCA(mouse_3w, npcs = 20)
# Here you should save the file? or do a Pearson correlation analysis
cor(mouse_3w)

mouse_3w <- RunUMAP(mouse_3w, dims = 1:20)
mouse_3w <- IntegrateLayers(object = mouse_3w,
                            method = RPCAIntegration,
                            orig.reduction = "pca",
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_3w <- JoinLayers(mouse_3w)

save(mouse_3w) # Try to save as data frame or something that is easier to open

# Save Seurat object
SaveH5Seurat(mouse_3w, overwrite = TRUE)

## Start here to directly load the Seurat object
# Load Seurat object
hfile <- Connect("Full Merge Project.h5Seurat")
hfile


## Filter for ventricular cardiomyocytes
mouse_vcm_3w <- subset(mouse_3w,
                       predicted.celltype.l2.score >= 0.7 &
                         mapping.score >= 0.7 &
                         predicted.celltype.l1.score >= 0.7 &
                         mapping.score >= 0.7 &
                         predicted.celltype.l2 == 'Ventricular Cardiomycoyte')
head(mouse_vcm_3w, 5) #control post: inspect the top 5 rows in the dataset
tail(mouse_vcm_3w, 5)

SaveH5Seurat(mouse_vcm_3w, "V:/Silje/snRNAseq/3weeks_sham_AB/Mouse_vcm_3w", overwrite = TRUE) # Save Seurat object?


## Find top genes
# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(mouse_vcm_3w) <- mouse_vcm_3w@meta.data$orig.ident

# Find differentially expressed genes
markers <- FindMarkers(mouse_3w_vcm, ident.1 = "3 weeks - AB", ident.2 = "3 weeks - SHAM")

markers <- markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))

write_xlsx(all_genes_vcm_3w, "V:/Silje/snRNAseq/3weeks_sham_AB/all_genes_vcm_3w.xlsx") # Save data frame as excel


top_genes_vcm_3w <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC)))

write_xlsx(top_genes_vcm_3w, "V:/Silje/snRNAseq/3weeks_sham_AB/top_genes_vcm_3w.xlsx") # Save data frame as excel

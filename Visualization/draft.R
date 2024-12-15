# Combine all time points into one seurat object
# Silje Wilhelmsen


library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)


    
mouse_vcm_6h <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_6h.Rds")
mouse_vcm_12h <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_12h.Rds")
mouse_vcm_1d <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_1d.Rds")
mouse_vcm_3d <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_3d.Rds")
mouse_vcm_1w <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_1w.Rds")
mouse_vcm_3w <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_3w.Rds")


mouse_vcm_all_time_points <- merge(mouse_vcm_6h, 
                                   y = c(mouse_vcm_12h, 
                                         mouse_vcm_1d, 
                                         mouse_vcm_3d, 
                                         mouse_vcm_1w, 
                                         mouse_vcm_3w), 
                                   project = "mouse_Full_Merge_Project", 
                                   merge.data = TRUE)

mouse_vcm_all_time_points <- NormalizeData(mouse_vcm_all_time_points, 
                          normalization.method = "LogNormalize", 
                          scale.factor = 10000)
mouse_vcm_all_time_points <- FindVariableFeatures(mouse_vcm_all_time_points)
mouse_vcm_all_time_points <- ScaleData(mouse_vcm_all_time_points)
mouse_vcm_all_time_points <- RunPCA(mouse_vcm_all_time_points, npcs = 20)
mouse_vcm_all_time_points <- RunUMAP(mouse_vcm_all_time_points, dims = 1:20)
mouse_vcm_all_time_points <- IntegrateLayers(object = mouse_vcm_all_time_points, 
                            method = RPCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.rpca",
                            verbose = FALSE)
mouse_vcm_all_time_points_join_layers <- JoinLayers(mouse_vcm_all_time_points)

SaveH5Seurat(mouse_vcm_all_time_points_join_layers, 
             filename = "mouse_all_time_points_join_layers.h5seurat", 
             overwrite = TRUE)

saveRDS(mouse_vcm_all_time_points_join_layers)
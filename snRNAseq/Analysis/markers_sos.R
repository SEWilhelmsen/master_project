# Markers 

library(presto)
library(dplyr)

data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_6h_vcm_with_stress_status.Rds")
data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_12h_vcm_with_stress_status.Rds") # Done
data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_1d_vcm_with_stress_status.Rds") # Done
data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_3d_vcm_with_stress_status.Rds") # Done
data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_1w_vcm_with_stress_status.Rds") # Done
data <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_3w_vcm_with_stress_status.Rds") # Done



markers <- presto:::wilcoxauc.Seurat(X = data, group_by = 'Group', assay = 'data', seurat_assay = 'RNA')
markers <- markers %>% 
  mutate(Timepoint = "6 Hours") # Change time point
write.xlsx(markers, file = "markers_6h_vcm_group.xlsx")

markers <- presto:::wilcoxauc.Seurat(X = data, group_by = 'Stress_Status', assay = 'data', seurat_assay = 'RNA')
markers <- markers %>% 
  mutate(Timepoint = "6 Hours") # Change time point
write.xlsx(markers, file = "markers_6h_vcm_stress.xlsx")


# Combine average log fold change for all time points into one data frame:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_combine_logFC_for_all_timepoints.R")


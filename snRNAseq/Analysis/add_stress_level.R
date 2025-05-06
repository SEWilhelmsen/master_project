# Filter for stress

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
# library(Azimuth)
# library(SeuratData)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(harmony)



# Load data
mouse_vcm_all_time_points <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points.Rds")

# # Prepare addition of column in meta data
# sham_samples <- c('Sample1-TP75B', 'Sample3-TP52A', 'TP3-1', 'TP3-4', 'TP4-1', 
#                   'TP5-1', 'TP9-1', 'TP9-5', 'TP11-5', 'TP17-2','TP13-1','TP17-3','TP21-3','TP22-1','TP22-2','TP23-2','TP23-4','TP26-5')
# orab_samples <- c('TP14-2', 'TP11-1A', 'Sample2-TP511B', 'Sample4-TP62A', 'TP10-2', 
#                   'TP9-2', 'TP52-2', 'TP58-4', 'TP2-4', 'TP51-3','TP14-1','TP15-3','TP19-2', 'TP19-4','TP22-3','TP23-1','TP24-5','TP26-4')
# 
# # Define the conditions for each timepoint
# timepoint_06_hours <- c('TP19-2', 'TP19-4','TP21-3','TP22-1','TP22-2','TP22-3')
# timepoint_12_hours <- c('TP23-1', 'TP23-2','TP23-4','TP24-5','TP26-4','TP26-5')
# timepoint_24_hours <- c('TP17-2', 'TP14-2','TP13-1','TP14-1','TP17-3','TP15-3')
# timepoint_1_week <- c('TP2-4', 'TP3-1', 'TP3-4', 'TP4-1', 'TP52-2', 'TP58-4')
# timepoint_3_days <- c('TP9-1', 'TP9-2', 'TP9-5', 'TP10-2', 'TP11-5', 'TP11-1A')
# timepoint_3_week <- c('TP5-1','Sample1-TP75B', 'Sample3-TP52A','Sample2-TP511B', 'Sample4-TP62A','TP51-3')
# 
# # Create the new 'Timepoint' column and assign values based on the conditions
# mouse_vcm_all_time_points@meta.data$Timepoint <- ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_24_hours, '1 Day',
#                                                         ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_1_week, '1 Week',
#                                                                ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_3_days, '3 Days',
#                                                                       ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_06_hours, '6 Hours',
#                                                                              ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_12_hours, '12 Hours',
#                                                                                     '3 Weeks')))))


# Prepare addition of column in meta data
sham_samples <- c('6 Hours - SHAM', '12 Hours - SHAM', '1 Day - SHAM', '3 Days - SHAM', '1 week - SHAM', '3 Weeks - SHAM')
orab_samples <- c('6 Hours - AB', '12 Hours - AB', '1 Day - AB', '3 Days - AB', '1 week - AB', '3 Weeks - AB')

# Define the conditions for each timepoint
timepoint_06_hours <- c('6 Hours - SHAM', '6 Hours - AB')
timepoint_12_hours <- c('12 Hours - SHAM', '12 Hours - AB')
timepoint_24_hours <- c('1 Day - SHAM', '1 Day - AB')
timepoint_1_week <- c('3 Days - SHAM', '3 Days - AB')
timepoint_3_days <- c('1 week - SHAM', '1 week - AB')
timepoint_3_week <- c('3 Weeks - SHAM', '3 Weeks - AB')

# Create the new 'Timepoint' column and assign values based on the conditions
mouse_vcm_all_time_points@meta.data$Timepoint <- ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_24_hours, '1 Day',
                                                        ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_1_week, '1 Week',
                                                               ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_3_days, '3 Days',
                                                                      ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_06_hours, '6 Hours',
                                                                             ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% timepoint_12_hours, '12 Hours',
                                                                                    '3 Weeks')))))



# Create a new column in meta data called Group
mouse_vcm_all_time_points@meta.data$Group <- ifelse(mouse_vcm_all_time_points@meta.data$orig.ident %in% sham_samples, 'SHAM', 'ORAB')

## Create the new group and timepoint 
mouse_vcm_all_time_points@meta.data$Group_Timepoint <- paste(mouse_vcm_all_time_points@meta.data$Group, mouse_vcm_all_time_points@meta.data$Timepoint, sep = "_")


# Verify outcomes
head(mouse_vcm_all_time_points@meta.data)
tail(mouse_vcm_all_time_points@meta.data)
print(unique(mouse_vcm_all_time_points@meta.data$Group))


# Genes used for custom definition of stress
genes_of_interest <- c("MYH7", "NPPA", "NPPB", "ANKRD1")

DefaultAssay(mouse_vcm_all_time_points) <- "RNA" # Which is the correct one? 
thresholds <- apply(FetchData(mouse_vcm_all_time_points, vars = genes_of_interest), 2, function(x) mean(x, na.rm = TRUE))

View(thresholds)

# If cardiomyocytes express 5x mean expression of all genes, they are characterized as stressed
selected_nuclei <- WhichCells(mouse_vcm_all_time_points, expression = MYH7 > 5* thresholds["MYH7"] |
                                NPPA > 5*  thresholds["NPPA"] |
                                NPPB >5*  thresholds["NPPB"] |
                                ANKRD1 > 5* thresholds["ANKRD1"])

meta_data <- mouse_vcm_all_time_points@meta.data

head(mouse_vcm_all_time_points@meta.data)

# Ensure the necessary columns exist
if (!all(c("predicted.celltype.l2", "Group") %in% colnames(meta_data))) {
  stop("Required columns not found in metadata")
}

# Create a new column for classification
meta_data$Stress_Status <- "Not classified"  # Default value

# Assign 'Stressed CM' if criteria are met
meta_data$Stress_Status[meta_data$predicted.celltype.l2 == "Ventricular Cardiomycoyte" & 
                          rownames(meta_data) %in% selected_nuclei & 
                          meta_data$Group == "ORAB"] <- "Stressed CM"

# Assign 'Non-CM' if the nucleus is not a Ventricular Cardiomyocyte
meta_data$Stress_Status[meta_data$predicted.celltype.l2 != "Ventricular Cardiomycoyte"] <- "Non-CM"

# Assign 'Not stressed CM' if it's a Ventricular Cardiomyocyte but not in the selected nuclei (Group AB)
meta_data$Stress_Status[meta_data$predicted.celltype.l2 == "Ventricular Cardiomycoyte" & 
                          !(rownames(meta_data) %in% selected_nuclei) & 
                          meta_data$Group == "ORAB"] <- "Not stressed CM"

# Assign 'SHAM - CM' if the nucleus belongs to Group SHAM
meta_data$Stress_Status[meta_data$predicted.celltype.l2 == "Ventricular Cardiomycoyte" & 
                          meta_data$Group == "SHAM"] <- "SHAM - CM"

# Assign updated metadata back to the Seurat object
mouse_vcm_all_time_points@meta.data <- meta_data
head(mouse_vcm_all_time_points@meta.data)
unique(mouse_vcm_all_time_points@meta.data$Stress_Status)

saveRDS(mouse_vcm_all_time_points, file = "mouse_vcm_all_time_points_with_stress_status.Rds")


# Subset for stressed data 
mouse_vcm_stressed <- subset(mouse_vcm_all_time_points, cells = selected_nuclei)
saveRDS(mouse_vcm_stressed, file = "mouse_vcm_stressed.Rds")


# Load data
mouse_vcm_stressed <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_stressed.Rds")


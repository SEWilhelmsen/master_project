# Visualize spots with expression of Nppa, Ankrd1, Myh7 and Nppb

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(spacexr)
library(tidyr)


# Load data and add columns
# data_object <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/data_object.Rds")

# ## 1 day ### 
# ### AB ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP39.4_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP39-4-1/outs"
prefix <- "TP39.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP39.4"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "1 Day"


# ## 1 day ### 
# ### SHAM ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP38.4_RCTD.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP38-4-2/outs"
prefix <- "TP38.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP38.4"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "1 Day"


# ## 3 days ### 
# ### ORAB ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP36.1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP36-1-1/outs"
prefix <- "TP36.1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP36.1"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Days"

RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP37.4_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP37-4-1/outs"
prefix <- "TP37.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP37.4"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Days"

RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP37.2_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP37-2/outs"
prefix <- "TP37.2_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP37.2"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Days"


# ## 3 days ### 
# ### SHAM ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP37.5_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP37-5-3/outs"
prefix <- "TP37.5_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP37.5"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Days"

RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP36.3_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP36-3/outs"
prefix <- "TP36.3_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP36.3"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Days"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP36.4_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP36-4/outs"
prefix <- "TP36.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP36.4"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Days"


## 1 Week ### 
### ORAB ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP52.4_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP52-4/outs"
prefix <- "TP52.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP52.4"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "1 Week"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP52.3_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP52-3/outs"
prefix <- "TP52.3_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP52.3"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "1 Week"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP46.4_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP46-4/outs"
prefix <- "TP46.4_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP46.4"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "1 Week"


## 1 Week ### 
### SHAM ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP46.2_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP46-2/outs"
prefix <- "TP46.2_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP46.2"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "1 Week"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP42.1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-1/outs"
prefix <- "TP42.1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP42.1"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "1 Week"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP42.2_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-2/outs"
prefix <- "TP42.2_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP42.2"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "1 Week"


## 3 Week ### 
### ORAB ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP32.3_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP32-3-3/outs"
prefix <- "TP32.3_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP32.3"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Weeks"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/A1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/A1/outs"
prefix <- "A1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "A1"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Weeks"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/D2_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D2/outs"
prefix <- "D2_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "D2"
object@meta.data$group <- "AB"
object@meta.data$Timepoint <- "3 Weeks"


## 3 Weeks ### 
### SHAM ###
RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/D1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D1/outs"
prefix <-"D1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "D1"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Weeks"

RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP31.1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP31-1/outs"
prefix <- "TP31.1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP31.1"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Weeks"


RCTD <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/tmp/TP48.1_RCTD_temp.rds")
localdir <- "//hypatia.uio.no/lh-med-imb-cardiodata/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP48-1-1/outs"
prefix <- "TP48.1_"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))
object@meta.data$sample <- "TP48.1"
object@meta.data$group <- "SHAM"
object@meta.data$Timepoint <- "3 Weeks"





# For plotting
xlim_all <- c(-2000, 11000)
ylim_all <- c(-300, 11000)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)

data_mat <- GetAssayData(object, layer = "data")   # rader = gener, kolonner = barcodes


# NPPA 
#-------------------------------------------------------------------
if (!("Nppa" %in% rownames(data_mat))) stop("Gen 'Nppa' finnes ikke i data-matrisen.")

# Retrieve Nppa expression per barcode 
nppa_vals <- as.numeric(data_mat["Nppa", ])
names(nppa_vals) <- colnames(data_mat)

# # Match values by barcodes 
object@meta.data$nppa_expr <- NA
matching <- intersect(rownames(object@meta.data), names(nppa_vals))
object@meta.data[matching, "nppa_expr"] <- nppa_vals[matching]


# Estimate threshold
# threshold_nppa <- mean(nppa_vals, na.rm = TRUE) * 5
# high_vec_nppa <- nppa_vals > threshold_nppa  


# Retrieve threshold of SHAM
threshold_sham_nppa <- thresholds_stress %>%
  filter(group == "SHAM", timepoint == "3 Weeks") %>%
  pull(threshold_nppa)
print(threshold_sham_nppa)
threshold_sham_nppa <- mean(threshold_sham_nppa, na.rm = TRUE)
print(threshold_sham_nppa)

high_vec_nppa <- nppa_vals > threshold_sham_nppa  


rctd_bcs <- colnames(RCTD@spatialRNA@counts)
high_for_rctd <- setNames(rep(FALSE, length(rctd_bcs)), rctd_bcs)
common <- intersect(names(high_vec_nppa), rctd_bcs)
high_for_rctd[common] <- high_vec_nppa[common]

# Investigate number of Nppa spots
cat("Number of high Nppa spots (in RCTD barcodes):", sum(high_for_rctd), "\n")

# Create coords_df for plotting
coords <- as.data.frame(RCTD@spatialRNA@coords)
coords_df <- data.frame(
  x = coords[,1],
  y = coords[,2],
  barcode = rownames(coords),
  high_nppa = high_for_rctd[rownames(coords)]
)
# Handle NAs as FALSE 
coords_df$high_nppa[is.na(coords_df$high_nppa)] <- FALSE

# Visualize
plot_nppa <- ggplot(coords_df, aes(x = x, y = y)) +
  geom_point(data = subset(coords_df, !high_nppa), color = "lightgrey", size = 1.2) +
  geom_point(data = subset(coords_df, high_nppa), color = "red", size = 1.6) +
  coord_equal() +
  xlim(xlim_all) + ylim(ylim_all) +
  theme_minimal() +
  labs(title = "Spots with expression of Nppa > 5x mean of SHAM")

plot_nppa

# ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "nppa.tiff"), 
#        plot = plot_nppa, width = 15, height = 15, dpi = 600)

ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "nppa_mean_sham.tiff"), 
       plot = plot_nppa, width = 15, height = 15, dpi = 600)



# Ankrd1 
#-------------------------------------------------------------------
if (!("Ankrd1" %in% rownames(data_mat))) stop("Gen 'Ankrd1' finnes ikke i data-matrisen.")

# Retrieve Ankrd1 expression per barcode 
ankrd1_vals <- as.numeric(data_mat["Ankrd1", ])
names(ankrd1_vals) <- colnames(data_mat)

# # Match values by barcodes 
object@meta.data$ankrd1_expr <- NA
matching <- intersect(rownames(object@meta.data), names(ankrd1_vals))
object@meta.data[matching, "ankrd1_expr"] <- ankrd1_vals[matching]


# Estimate threshold
# threshold_ankrd1 <- mean(ankrd1_vals, na.rm = TRUE) * 5
# high_vec_ankrd1 <- ankrd1_vals > threshold_ankrd1  


# Retrieve threshold of SHAM
threshold_sham_ankrd1 <- thresholds_stress %>%
  filter(group == "SHAM", timepoint == "3 Weeks") %>%
  pull(threshold_ankrd1)
print(threshold_sham_ankrd1)
threshold_sham_ankrd1 <- mean(threshold_sham_ankrd1, na.rm = TRUE)
print(threshold_sham_ankrd1)

high_vec_ankrd1 <- ankrd1_vals > threshold_sham_ankrd1  

rctd_bcs <- colnames(RCTD@spatialRNA@counts)
high_for_rctd <- setNames(rep(FALSE, length(rctd_bcs)), rctd_bcs)
common <- intersect(names(high_vec_ankrd1), rctd_bcs)
high_for_rctd[common] <- high_vec_ankrd1[common]

# Investigate number of Ankrd1 spots
cat("Number of high Ankrd1 spots (in RCTD barcodes):", sum(high_for_rctd), "\n")

# Create coords_df for plotting
coords <- as.data.frame(RCTD@spatialRNA@coords)
coords_df <- data.frame(
  x = coords[,1],
  y = coords[,2],
  barcode = rownames(coords),
  high_ankrd1 = high_for_rctd[rownames(coords)]
)
# Handle NAs as FALSE 
coords_df$high_ankrd1[is.na(coords_df$high_ankrd1)] <- FALSE

# Visualize
plot_ankrd1 <- ggplot(coords_df, aes(x = x, y = y)) +
  geom_point(data = subset(coords_df, !high_ankrd1), color = "lightgrey", size = 1.2) +
  geom_point(data = subset(coords_df, high_ankrd1), color = "red", size = 1.6) +
  coord_equal() +
  xlim(xlim_all) + ylim(ylim_all) +
  theme_minimal() +
  labs(title = "Spots with expression of Ankrd1 > 5x mean of SHAM")
  
plot_ankrd1

# ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "ankrd1.tiff"), 
#        plot = plot_ankrd1, width = 15, height = 15, dpi = 600)


ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "ankrd1_mean_sham.tiff"), 
       plot = plot_ankrd1, width = 15, height = 15, dpi = 600)



# Myh7 
#-------------------------------------------------------------------
if (!("Myh7" %in% rownames(data_mat))) stop("Gen 'Myh7' finnes ikke i data-matrisen.")

# Retrieve Myh7 expression per barcode 
myh7_vals <- as.numeric(data_mat["Myh7", ])
names(myh7_vals) <- colnames(data_mat)

# # Match values by barcodes 
object@meta.data$myh7_expr <- NA
matching <- intersect(rownames(object@meta.data), names(myh7_vals))
object@meta.data[matching, "myh7_expr"] <- myh7_vals[matching]


# # Estimate threshold
# threshold_myh7 <- mean(myh7_vals, na.rm = TRUE) * 5
# high_vec_myh7 <- myh7_vals > threshold_myh7  

# Retrieve threshold of SHAM
threshold_sham_myh7 <- thresholds_stress %>%
  filter(group == "SHAM", timepoint == "3 Weeks") %>% # Change time point
  pull(threshold_myh7)
print(threshold_sham_myh7)
threshold_sham_myh7 <- mean(threshold_sham_myh7, na.rm = TRUE)
print(threshold_sham_myh7)

high_vec_myh7 <- myh7_vals > threshold_sham_myh7  


rctd_bcs <- colnames(RCTD@spatialRNA@counts)
high_for_rctd <- setNames(rep(FALSE, length(rctd_bcs)), rctd_bcs)
common <- intersect(names(high_vec_myh7), rctd_bcs)
high_for_rctd[common] <- high_vec_myh7[common]

# Investigate number of Myh7 spots
cat("Number of high Myh7 spots (in RCTD barcodes):", sum(high_for_rctd), "\n")

# Create coords_df for plotting
coords <- as.data.frame(RCTD@spatialRNA@coords)
coords_df <- data.frame(
  x = coords[,1],
  y = coords[,2],
  barcode = rownames(coords),
  high_myh7 = high_for_rctd[rownames(coords)]
)
# Handle NAs as FALSE 
coords_df$high_myh7[is.na(coords_df$high_myh7)] <- FALSE

# Visualize
plot_myh7 <- ggplot(coords_df, aes(x = x, y = y)) +
  geom_point(data = subset(coords_df, !high_myh7), color = "lightgrey", size = 1.2) +
  geom_point(data = subset(coords_df, high_myh7), color = "red", size = 1.6) +
  coord_equal() +
  xlim(xlim_all) + ylim(ylim_all) +
  theme_minimal() +
  labs(title = "Spots with expression of Myh7 > 5x mean of SHAM")

plot_myh7

# ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "myh7.tiff"), 
#        plot = plot_myh7, width = 15, height = 15, dpi = 600)

ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "myh7_mean_sham.tiff"), 
       plot = plot_myh7, width = 15, height = 15, dpi = 600)


# Nppb 
#-------------------------------------------------------------------
if (!("Nppb" %in% rownames(data_mat))) stop("Gen 'Nppb' finnes ikke i data-matrisen.")

# Retrieve Nppb expression per barcode 
nppb_vals <- as.numeric(data_mat["Nppb", ])
names(nppb_vals) <- colnames(data_mat)

# # Match values by barcodes 
object@meta.data$nppb_expr <- NA
matching <- intersect(rownames(object@meta.data), names(nppb_vals))
object@meta.data[matching, "nppb_expr"] <- nppb_vals[matching]


# # Estimate threshold
# threshold_nppb <- mean(nppb_vals, na.rm = TRUE) * 5
# high_vec_nppb <- nppb_vals > threshold_nppb  

# Retrieve threshold of SHAM
threshold_sham_nppb <- thresholds_stress %>%
  filter(group == "SHAM", timepoint == "3 Weeks") %>%
  pull(threshold_nppb)
print(threshold_sham_nppb)
threshold_sham_nppb <- mean(threshold_sham_nppb, na.rm = TRUE)
print(threshold_sham_nppb)

high_vec_nppb <- nppb_vals > threshold_sham_nppb  


rctd_bcs <- colnames(RCTD@spatialRNA@counts)
high_for_rctd <- setNames(rep(FALSE, length(rctd_bcs)), rctd_bcs)
common <- intersect(names(high_vec_nppb), rctd_bcs)
high_for_rctd[common] <- high_vec_nppb[common]

# Investigate number of nppb spots
cat("Number of high nppb spots (in RCTD barcodes):", sum(high_for_rctd), "\n")

# Create coords_df for plotting
coords <- as.data.frame(RCTD@spatialRNA@coords)
coords_df <- data.frame(
  x = coords[,1],
  y = coords[,2],
  barcode = rownames(coords),
  high_nppb = high_for_rctd[rownames(coords)]
)
# Handle NAs as FALSE 
coords_df$high_nppb[is.na(coords_df$high_nppb)] <- FALSE

# Visualize
plot_nppb <- ggplot(coords_df, aes(x = x, y = y)) +
  geom_point(data = subset(coords_df, !high_nppb), color = "lightgrey", size = 1.2) +
  geom_point(data = subset(coords_df, high_nppb), color = "red", size = 1.6) +
  coord_equal() +
  xlim(xlim_all) + ylim(ylim_all) +
  theme_minimal() +
  labs(title = "Spots with expression of Nppb > 5x mean of SHAM")

plot_nppb

# ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "nppb.tiff"), 
#        plot = plot_nppb, width = 15, height = 15, dpi = 600)


ggsave(filename = paste0("C:/Users/siljeew/Master_project/spatial_transcriptomic/Plots/", prefix, "nppb_mean_sham.tiff"), 
       plot = plot_nppb, width = 15, height = 15, dpi = 600)



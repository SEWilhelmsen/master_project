# Create spatials plots for library of cell types

library(spacexr)
library(Seurat)
library(ggplot2) #for saving result plots
## https://satijalab.org/seurat/articles/visiumhd_commands_intro#visium-hd-support-in-seurat
##https://www.10xgenomics.com/analysis-guides/integrating-10x-visium-and-chromium-data-with-r
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


load("C:/Users/Labuser/master_project/SCreference_l2_ALLNUCLEI_10112025.RData")

# ## 1 day ### 
# ### ORAB ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP39-4-1/outs"
prefix <- "TP39.4_"
# ## 1 day ### 
# ### SHAM ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP38-4-2/outs"
prefix <- "TP38.4_"
# 
# ## 3 days ### 
# ### ORAB ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP36-1-1/outs"
prefix <- "TP36.1_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP37-4-1/outs"
prefix <- "TP37.4_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP37-2/outs"
prefix <- "TP37.2_"

# ## 3 days ### 
# ### SHAM ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP37-5-3/outs"
prefix <- "TP37.5_"
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP36-3/outs"
prefix <- "TP36.3_"
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP36-4/outs"
prefix <- "TP36.4_"


## 1 Week ### 
### ORAB ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP52-4/outs"
prefix <- "TP52.4_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP52-3/outs"
prefix <- "TP52.3_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D2-TP46-4/outs"
prefix <- "TP46.4_"

## 1 Week ### 
### SHAM ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A2-TP46-2/outs"
prefix <- "TP46.2_"
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-1/outs"
prefix <- "TP42.1_"
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP42-2/outs"
prefix <- "TP42.2_"

## 3 Week ### 
### ORAB ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/D1-TP32-3-3/outs"
prefix <- "TP32.3_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/A1/outs"
prefix <- "A1_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D2/outs"
prefix <- "D2_"

## 3 Weeks ### 
### SHAM ###
localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0901/GCF0901/D1/outs"
prefix <-"D1_"

# localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP31-1/outs"
localdir <- "C:/Users/Labuser/master_project/Spatials/tmp/A1-TP31-1/outs"
prefix <- "TP31.1_"

localdir <- "D:/CardioTarget/Spatial Transcriptomics/GCF-0902/SpaceRanger/A1-TP48-1-1/outs"
prefix <- "TP48.1_"

# Load the data
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# SCTransform also supported, but we use LogNormalize here
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)




DefaultAssay(object) <- "Spatial.016um"
counts_hd <- object[["Spatial.016um"]]$counts
object_cells_hd <- colnames(object[["Spatial.016um"]])
coords <- GetTissueCoordinates(object)[object_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

rownames(SCreference@counts) <- toupper(rownames(SCreference@counts))
rownames(query@counts) <- toupper(rownames(query@counts))

RCTD <- create.RCTD(query, SCreference, max_cores = 6)
RCTD <- run.RCTD(RCTD, doublet_mode = "full") # this command takes ~15 mins to run

barcodes <- colnames(RCTD@spatialRNA@counts) # list of spatial barcodes
weights <- RCTD@results$weights # Weights for each cell type per barcode

norm_weights <- normalize_weights(weights)
cell_type_names<-colnames(norm_weights) # List of cell types
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df); rownames(subset_df) <- NULL


plot_dir <- "C:/Users/Labuser/master_project/Spatials/Plots"

# Sjekk at output-mappen finnes
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Sjekk at nødvendige objekter finnes
stopifnot(exists("RCTD"), exists("norm_weights"), exists("cell_type_names"), exists("barcodes"))

# Antall celletyper
n_types <- length(cell_type_names)  # eller n_types <- 19 hvis du vet det

# Hent coords — vi trenger den bare om du vil beregne limits automatisk
if (!is.null(RCTD@spatialRNA@coords)) {
  coords <- RCTD@spatialRNA@coords
} else if (!is.null(RCTD@spatialRNA@coordinates)) {
  coords <- RCTD@spatialRNA@coordinates
} else {
  try({
    coords <- GetTissueCoordinates(object)[colnames(RCTD@spatialRNA@counts), 1:2]
  }, silent = TRUE)
}
if (is.null(coords)) {
  warning("Finner ikke coords automatisk — du bruker hardkodede xlim_all / ylim_all.")
} else {
  coords <- as.matrix(coords)
  # hvis ønskelig, kan du beregne xlim_all / ylim_all her (men du har hardkodet dem)
}

# Hardkodede globale limits (bruk disse for å gjenbruke samme akser på andre samples)
xlim_all <- c(-2000, 11000)
ylim_all <- c(-300, 11000)

# Preview plot for rotation
i <- 19
plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[,cell_type_names[i]], title =cell_type_names[i], size=1)

# # Save RCTD
# data_dir <- "C:/Users/Labuser/master_project/Spatials/tmp/"
# saveRDS(RCTD, file = paste0(data_dir, prefix, "RCTD.rds"))

# Rotate slides if necessary
# open: "C:/Users/Labuser/master_project/Spatials/Analysis/rotation_of_slide.R"


# # Test om plot_puck_continuous returnerer en ggplot
# test_plot <- try(plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[, cell_type_names[1]], title = cell_type_names[1], size = 1), silent = TRUE)
# ggplot_mode <- inherits(test_plot, "ggplot")
# 
# if (ggplot_mode) {
#   library(ggplot2)
#   for (i in seq_len(n_types)) {
#     cellname <- cell_type_names[i]
#     p <- plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[, cellname], title = cellname, size = 1)
#     p_fixed <- p + coord_fixed(xlim = xlim_all, ylim = ylim_all, expand = FALSE)
#     fname <- file.path(plot_dir, paste0(prefix, "distribution_", i, "_", cellname, ".tiff"))
#     ggsave(filename = fname, plot = p_fixed, device = "tiff", width = 10, height = 8, units = "in", dpi = 600)
#   }
# } else {
#   for (i in seq_len(n_types)) {
#     cellname <- cell_type_names[i]
#     fname <- file.path(plot_dir, paste0(prefix, "distribution_", i, "_", cellname, ".tiff"))
#     tiff(filename = fname, width = 10, height = 8, units = "in", res = 600)
#     op <- par(pty = "s", mar = c(3,3,3,3))
#     # Prøv å sende limits — hvis funksjonen ignorerer dem, fallback til å bare kalle funksjonen
#     try({
#       plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[, cellname],
#                            title = cellname, size = 1, xlim = xlim_all, ylim = ylim_all, asp = 1)
#     }, silent = TRUE)
#     # Kall uansett for å være sikker
#     tryCatch({
#       plot_puck_continuous(RCTD@spatialRNA, barcodes, norm_weights[, cellname], title = cellname, size = 1)
#     }, error = function(e) {
#       message("Feil ved plotting for ", cellname, ": ", e$message)
#     })
#     par(op)
#     dev.off()
#   }
# }


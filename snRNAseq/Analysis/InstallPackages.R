# Install packages 
options(timeout = 250) 

library(SoupX)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(cowplot)

# Install package: signac
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("stuart-lab/signac", ref = "develop")
library(Signac)


# Install package: GenomeInfoDb
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomeInfoDb")
library(GenomeInfoDb)

# Install package: EnsDb.Mmusculus.v79
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Mmusculus.v79")
library(EnsDb.Mmusculus.v79)

# Install package: biovizBase
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biovizBase")
library(biovizBase)

# Install package: Seurat
BiocManager::install("Seurat")
BiocManager::install(version = "3.18")
library(Seurat)

# Install package: TFBSTools
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TFBSTools", type = "source", force = TRUE)
library(TFBSTools)

# Install package: Azimuth
# First, install dependencies and check for updates
update.packages(oldPkgs = c("withr", "rlang"))

# Install azimuth
install.packages("remotes")
remotes::install_github("satijalab/azimuth")
library(Azimuth)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github('satijalab/azimuth', ref = 'master')

# Install package: SeuratData
devtools::install_github('satijalab/seurat-data')
library(SeuratData)

# Install reference data from SeuratData
install.packages("heartref.SeuratData", repos = "http://seurat.nygenome.org", type = "source")


# Install package: Spacexr
install.packages("devtools")
options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)

# Install package: Matrix
install.packages("Matrix")
install.packages("Matrix", repos = "http://R-Forge.R-project.org")
library(Matrix)

# Install package: zellkonverter
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("zellkonverter")
library(zellkonverter)

# Install package:SingleCellExperiment
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)

# Install package:reticulate
install.packages("reticulate")
library(reticulate)

# Install package: harmony
install.packages("harmony")
library(harmony)

# Install package: SeuratDisk
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

# Install package: somepackage
install.packages("remotes")
remotes::install_github("kirillseva/somepackage", force = TRUE)
library(somepackage)

# Install package: clusterProfiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)

# Install package: org.Hs.eg.db
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# Install package: pathview
# install from BioConductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("datapplab/pathview")

# Install package: dplyr
install.packages("dplyr")
library(dplyr)

# Install package: ggplot2
install.packages("ggplot2")
library(ggplot2)

# Install EnhancedVolcano
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)

# Plot EnhancedVolcano
EnhancedVolcano(
  markers,
  lab = markers$gene,
  x = 'Log2 Fold Change',
  y = '-Log10 Adjusted P-value',
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano Plot",
  subtitle = "Differential Expression Analysis"
)
head(markers)

# Install package: chromVAR
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("chromVAR")


# Install package: JASPAR2020
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("JASPAR2020")


# Install package: motifmatchr 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("motifmatchr")


# Install package: BSgenome.Mmusculus.UCSC.mm10
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")


# Update packages 
update.packages()


# Check library path
.libPaths()


# FC%r tabell
table(Prediction_TP511@meta.data$predicted.celltype.l2)

# Change size allowed to process 
options(future.globals.maxSize = 4 * 1024^3) # This sets the maximum size to 4 GB. 


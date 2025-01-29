## Install packages for ATAC

library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
#library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(Seurat)

if (!"BiocManager" %in% rownames(installed.packages()))
  install.packages("BiocManager")
BiocManager::install("motifTestR")

library(motifTestR)
library(extraChIPs)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(parallel)
library(ggplot2)
library(patchwork)
library(universalmotif)


library(dplyr)
library(patchwork)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(Seurat)
        
# TP172
counts_tp172 <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP17-2/outs/filtered_feature_bc_matrix.h5",
                           use.names = TRUE, unique.features = TRUE)
fragments_tp172 <- "W:/CardioTarget/Multiome seq data/TP17-2/outs/atac_fragments.tsv.gz"
TP172 <- CreateSeuratObject(
  counts = counts_tp172$Peaks,
  assay = "ATAC",
  project = "1 Day - SHAM"
)

# Add Chromatin Assay
TP172[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp172$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp172,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)

# TP131
counts_tp131 <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP13-1/outs/filtered_feature_bc_matrix.h5",
                           use.names = TRUE, unique.features = TRUE)
fragments_tp131 <- "W:/CardioTarget/Multiome seq data/TP13-1/outs/atac_fragments.tsv.gz"
TP131 <- CreateSeuratObject(
  counts = counts_tp131$Peaks,
  assay = "ATAC",
  project = "1 Day - SHAM"
)


# Add Chromatin Assay
TP131[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp131$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp131,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)

# TP173
counts_tp173 <- Read10X_h5("W:/CardioTarget/Multiome seq data/TP17-3/outs/filtered_feature_bc_matrix.h5",
                           use.names = TRUE, unique.features = TRUE)
fragments_tp173 <- "W:/CardioTarget/Multiome seq data/TP17-3/outs/atac_fragments.tsv.gz"
TP173 <- CreateSeuratObject(
  counts = counts_tp173$Peaks,
  assay = "ATAC",
  project = "1 Day - SHAM"
)


# Add Chromatin Assay
TP173[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp173$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp173,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)


# TP142
counts_tp142 <- Read10X_h5("D:/MultiomeSeq/TP14-2/outs/filtered_feature_bc_matrix.h5", 
                           use.names = TRUE, unique.features = TRUE)
fragments_142 <- "D:/MultiomeSeq/TP14-2/outs/atac_fragments.tsv.gz"
TP142 <- CreateSeuratObject(
  counts = counts_tp142$Peaks,
  assay = "ATAC",
  project = "1 Day - AB"
)

# Add Chromatin Assay
TP142[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp142$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp142,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)

# TP141
counts_tp141 <- Read10X_h5("D:/MultiomeSeq/TP14-1/outs/filtered_feature_bc_matrix.h5",
                           use.names = TRUE, unique.features = TRUE)
fragments_141 <- "D:/MultiomeSeq/TP14-1/outs/atac_fragments.tsv.gz"
TP141 <- CreateSeuratObject(
  counts = counts_tp141$Peaks,
  assay = "ATAC",
  project = "1 Day - AB"
)

# Add Chromatin Assay
TP141[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp141$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp141,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)

# TP153
counts_tp153 <- Read10X_h5("D:/MultiomeSeq/TP15-3/outs/filtered_feature_bc_matrix.h5",
                           use.names = TRUE, unique.features = TRUE)
fragments_153 <- "D:/MultiomeSeq/TP15-3/outs/atac_fragments.tsv.gz"
TP153 <- CreateSeuratObject(
  counts = counts_tp153$Peaks,
  assay = "ATAC",
  project = "1 Day - AB"
)

# Add Chromatin Assay
TP153[['ATAC']] <- CreateChromatinAssay(
  counts = counts_tp153$Peaks,
  sep = c(":", "-"),
  fragments = fragments_tp153,
  annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79),
  genome = 'mm10'
)

# Optionally Merge or Integrate Seurat Objects
# seurat_combined <- merge(TP131, y = list(TP173, TP142))



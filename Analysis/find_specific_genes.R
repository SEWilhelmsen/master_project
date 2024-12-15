### Find only specific genes 
# Input: data frame of genes (rows) and timepoints (columns) based on avg_log2FC
## Output: data frame of only specific genes
# Silje Wilhelmsen


library(dplyr)
library(readr)

# Create a list of specific genes  
### Glycolysis 
# hexokinase or glucokinase: HK1, HK2, GCK,  
# phosphofructokinase: PFK1, PFKM, 
# pyruvate kinase: PK, 
#gene_list_glycolysis <- list("HK1", "HK2", "GCK", "PFK1", "PFKM", "PK")

### Binding from glycolysis to citric acid cycle
# pyruvate dehydrogenase complex: PDH, PDHA1, DLAT, 
#gene_list_pdh <- list("PDH", "PDHA1", "DLAT")

### tricarboxylic cycle 
# citrate synthase: CS, 
# oxoglutarate (a-ketoglutarate) dehydrogenase complex: OGDH, OGDC, DLST, DLD, 
# isocitrate dehydrogenase complex: IDH1, IDH2, IDH3A, IDH3G, IDH3B
#gene_list_tca <- list("CS", "OGDH", "OGDC", "DLST", "DLD", "IDH1", "IDH2", "IDH3A", "IDH3G", "IDH3B")


# Convert gene_list to uppercase to ensure uniform formatting
gene_list <- list(
  glycolysis = toupper(c("HK1", "HK2", "GCK", "PFK1", "PFKM", "PK")),
  pdh = toupper(c("PDH", "PDHA1", "PDHA2", "DLAT", "DLD", "PDK1", "PDK2",  "PDHX")),
  tca = toupper(c("CS", "OGDH", "OGDC", "DLST", "DLD", "IDH1", "IDH2", "IDH3A", "IDH3G", "IDH3B")),
  markers = toupper(c("ACTC1", "TPM", "RYR2",	"MYH6", "ATP2A2", "NPPA", "TNNC1", "ACTA1")),
  chylomicron_clearance = toupper(c("APOC1", "APOC2", "APOC3", "APOE", "LIPC", "GPIHBP")),
  uncoupling_proteins = toupper(c("SLC25A4", "SLC25A5", "UCP1", "UCP2", "UCP3")),
  fatty_acid_transport = toupper(c("FABP3", "FABPH", "CD36", "FATP4", "CACP", "MCAT")),
  b_oxidation = toupper(c("CPT1B", "CPT1A", "ACADVL", "CPT2"))
)

specific_genes <- unique(unlist(gene_list))
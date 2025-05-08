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

### Binding from glycolysis to citric acid cycle
# pyruvate dehydrogenase complex: PDH, PDHA1, DLAT, DLD, PDHX

### tricarboxylic cycle 
# citrate synthase: CS, 
# oxoglutarate (a-ketoglutarate) dehydrogenase complex: OGDH, OGDC, DLST, DLD? 
# isocitrate dehydrogenase complex: IDH1, IDH2, IDH3A, IDH3G, IDH3B

### Ketone catabolism
#   ketone_catabolism: "ACAT1", "OXCT1", "OXCT2", "OXCT2A", "OXCT2B", "SLC16A1" for transport


# Convert gene_list to uppercase to ensure uniform formatting
gene_list <- list(
  glycolysis = toupper(c("HK1", "HK2", "GCK", "PFKM")),
  pdh = toupper(c("PDHA1", "DLD", "DLAT",  "PDHX")),
  tca = toupper(c("IDH1", "IDH2", "IDH3B")),
  cm_markers = toupper(c("ACTC1", "TPM", "RYR2",	"MYH6", "ATP2A2", "NPPA", "TNNC1", "ACTA1")),
  chylomicron_clearance = toupper(c("APOC1", "APOC2", "APOC3", "APOE", "LIPC", "GPIHBP")),
  uncoupling_proteins = toupper(c("UCP1", "UCP2", "UCP3")),
  glucose_transmembrane_transport = toupper(c("SLC2A4")),
  fatty_acid_transport = toupper(c("CD36")),
  b_oxidation = toupper(c("CPT1A", "CPT1B", "CPT2")),
  regulation_of_ketone_metabolic_process = toupper(c("AVP", "SLC45A3", "PLIN5")),
  hydroxybutyrate_dehydrogenase_activity = toupper(c("BDH1")),
  ketone_catabolism = toupper(c("OXCT1", "OXCT2A", "OXCT2B")),
  protein_transmembrane_transport = toupper(c("")),
  cardiac_muscle_contraction = toupper(c("ATP2A2", "TNNI3","ATP1A1", "ATP1A2", "ATP1B1", "ATP2A1")),
  negative_regulation_of_cardiac_muscle_contraction = toupper(c("ADCY10", "BIN1", "PDE5A", "SRI", "ZC3H12A", "PIK3CG","MIR3")),
  positive_regulation_of_cardiac_muscle_contraction = toupper(c("ADRA1A", "HSP90AA1", "KCNQ1", "RGS2")),
  regulation_of_cardiac_calcium_signaling = toupper(c("ATP1A2", "ATP1B1", "ATP2A2", "CALM1", "CALM2", "CALM3", "PLN", "TNNI3")),
  stress_markers = toupper(c("NPPA", "NPPB", "MYH7", "ANKRD1")),
  negative_regulation_of_creb_transcription_factor_activity = toupper(c("DDIT3", "SIK1", "ADGRG3", "EIF2AK4")),
  positive_regulation_of_creb_transcription_factor_activity = toupper(c("ADCY8", "CD200", "RELN", "LRP8", "RPS6KA5", "CRTC1", "CAMK1D", "CRTC3", "CRTC2")),
  negative_regulation_of_nf_kappab_transcription_factor_activity = toupper(c("NFKBIA", "NFKBIB", "IRAK1", "IRAK2", "IRAK3")),
  nf_kappab_complex = toupper(c("NFKB1", "REL", "RELA", "RELB")),
  peroxisome_proliferator_activated_receptor_binding = toupper(c("DUT", "HMGA1", "PRMT2", "MDM2", "NFATC4", "MED1", "TACC1", "NR0B2", "ASXL2", "ASXL3", "ASXL1")),
  spec_genes = toupper(c("MYH7", "NPPA", "NPPB", "ANKRD1", "HK1", "HK2", "GCK", "PFK1", "PFKM", "PK"))
)

specific_genes <- unique(unlist(gene_list))
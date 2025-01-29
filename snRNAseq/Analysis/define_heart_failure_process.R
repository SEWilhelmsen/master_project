# Contractility in CM


# Load libraries
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(readr)

# Define gene ontology process
go_id_process <- list(
  cardiac_muscle_contraction = "GO:0060048",
  calcium_signaling_pathway = "GO:0019722",
  cardiac_muscle_tissue_growth = "GO:0055017"
)


# Function to retrieve genes for a GO term
get_genes_for_go_process <- function(go_term) {
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = go_term, 
                                      keytype = "GOALL", 
                                      columns = "ENTREZID")[, "ENTREZID"]
  
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = entrez_ids, 
                                        keytype = "ENTREZID", 
                                        columns = "SYMBOL")[, "SYMBOL"]
  
  # Normalize gene symbols
  gene_symbols_upper <- toupper(gene_symbols)
  
  # Return unique gene symbols related to the GO term
  unique(gene_symbols_upper)
}

# Function to retrieve genes organized by process
get_genes_of_interest <- function(go_terms, define_go_process, get_genes_for_go_process) {
  # Initialize a named list to store genes for each go process
  genes_list <- list()
  
  # For each go process, retrieve the genes and store into genes_list
  for (process in go_terms) {
    go_term <- define_go_process[[process]]
    if (is.null(go_term)) {
      warning(paste("GO term not found for the process:", process))
    } else {
      genes <- get_genes_for_go_process(go_term)
      genes_list[[process]] <- genes # Store list by process name
    }
  }
  
  print(genes_list)  # Debug print to check the expected structure
  
  genes_list  # Return the full list organized by process
}


# Convert gene_list to uppercase to ensure uniform formatting
gene_list <- list(
  cardiac_muscle_contraction = toupper(c("ADORA1", "ADRA1B", "ADRA1A", "GRK2", "BIN1", 
  "ANK2", "ATP1A1", "ATP1A2", "ATP1B1", "ATP2A1", "ATP2A2", "CACNA1C", "CACNA1D", "CACNA2D1", 
  "CACNB2", "CALM1", "CALM2", "CALM3", "CAMK2D", "CASQ2", "CAV1", "CAV3", "CHGA", "CLIC2", "NKX2-5", 
  "DLG1", "DMD", "DSC2", "DSG2", "DSP", "FGF12", "FGF13", "FKBP1B", "FLNA", "MTOR", "GAA", "GATA4", 
  "GJA1", "GJA5", "GSN", "GSTM2", "HRC", "HSP90AA1", "JUP", "KCNA5", "KCND3", "KCNE1", "KCNH2", "KCNJ2", 
  "KCNJ3", "KCNJ5", "KCNJ8", "KCNN2", "KCNQ1", "SMAD5", "SMAD7", "MYBPC3", "MYH6", "MYH7", "MYL2", "MYL3", 
  "MYL4", "NOS1", "NPPA", "P2RX4", "PDE4B", "PDE4D", "PIK3CA", "PIK3CG", "PKP2", "PLN", "PRKACA", "MAP2K3", 
  "MAP2K6", "RGS2", "RYR2", "SCN1A", "SCN1B", "SCN2B", "SCN4B", "SCN5A", "SCN10A", "SGCD", "SLC8A1", "SLC9A1", 
  "SNTA1", "SRI", "STC1", "TAFAZZIN", "TNNC1", "TNNI1", "TNNI2", "TNNI3", "TNNT2", "TPM1", "TTN", "SUMO1", "UCN", 
  "VEGFB", "CXCR4", "CSRP3", "TCAP", "PDE5A", "CACNA1G", "GSTO1", "NUP155", "NOS1AP", "KCNE2", "KCNE3", "HCN4", 
  "GJC1", "ABCC9", "AKAP9", "PPP1R13L", "GPD1L", "NEDD4L", "KCNE5", "KCNE4", "BMP10", "RANGRF", "CTNNA3", "EHD3", 
  "TNNI3K", "ASB3", "TRPM4", "TMEM38B", "SCN3B", "ADCY10", "ACE2", "TMEM38A", "ZC3H12A", "MYLK2", "RNF207", "MIR1-1", 
  "MIR133A1", "MIR200C", "MIR30E", "MIR328", "MIR448")),
  calcium_signaling_pathway = toupper(c()),
  cardiac_muscle_tissue_growth = toupper(c()),
  stress_markers = toupper(c("NPPA", "NPPB", "MYH7"))
)

specific_genes <- unique(unlist(gene_list))


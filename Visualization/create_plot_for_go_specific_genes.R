## Create plot of expression of genes in a specific list over time
# Input: mouse_vcm_all_genes_avg_log2fc
# Output: plot of gene expression over time
# Silje Wilhelmsen

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source("find_specific_genes.R")

# Convert gene_list to uppercase to ensure uniform formatting
gene_list <- list(
  glycolysis = toupper(c("HK1", "HK2", "GCK", "PFK1", "PFKM", "PK")),
  pdh = toupper(c("PDH", "PDHA1", "DLAT")),
  tca = toupper(c("CS", "OGDH", "OGDC", "DLST", "DLD", "IDH1", "IDH2", "IDH3A", "IDH3G", "IDH3B")),
  markers = toupper(c("ACTC1", "TPM", "RYR2",	"MYH6", "ATP2A2", "NPPA", "TNNC1", "ACTA1"))
)

#specific_genes <- unique(unlist(gene_list))

## Import data
mouse_vcm_all_genes_avg_log2fc <- read_csv("C:/Users/siljeew/snRNAseq/Data/mouse_vcm_all_genes_avg_log2fc.csv")


# Ensure gene names are uppercase
mouse_vcm_all_genes_avg_log2fc$gene <- toupper(mouse_vcm_all_genes_avg_log2fc$gene)

# Check the structure of the data
print(colnames(mouse_vcm_all_genes_avg_log2fc))
print(head(mouse_vcm_all_genes_avg_log2fc))


# Subset the data to include only genes in go_process_of_interest_genes
mouse_vcm_avg_log2fc_genes_of_interest <- mouse_vcm_all_genes_avg_log2fc %>%
  filter(gene %in% go_process_of_interest_genes)

# View the subsetted data
head(mouse_vcm_avg_log2fc_genes_of_interest)
View(mouse_vcm_avg_log2fc_genes_of_interest)  

# # Subset the data to include only the glycolysis genes
# glycolysis_genes <- gene_list$glycolysis

# Ensure gene names are uppercase in the dataset
mouse_vcm_avg_log2fc_genes_of_interest$gene <- toupper(mouse_vcm_avg_log2fc_genes_of_interest$gene)

# Filter for glycolysis genes
mouse_vcm_go_process_genes_avg_log2fc <- mouse_vcm_avg_log2fc_genes_of_interest %>%
  filter(gene %in% go_process_of_interest_genes)

# Reshape data from wide to long format
mouse_vcm_go_process_genes_avg_log2fc_long <- mouse_vcm_go_process_genes_avg_log2fc %>%
  pivot_longer(cols = -gene, names_to = "time_point", values_to = "avg_log2fc")

# Set the correct order of time points.
mouse_vcm_go_process_genes_avg_log2fc_long$time_point <- factor(mouse_vcm_go_process_genes_avg_log2fc_long$time_point,
                                                                levels = c("avg_log2fc_3d",
                                                                           "avg_log2fc_1w",
                                                                           "avg_log2fc_3w"))

# Print to confirm
print(head(mouse_vcm_go_process_genes_avg_log2fc_long))
nrow(mouse_vcm_go_process_genes_avg_log2fc_long)

# Create the ggplot for glycolysis genes
go_genes_plot <- ggplot(mouse_vcm_go_process_genes_avg_log2fc_long, aes(x = time_point, y = avg_log2fc, group = gene, color = gene)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Set1") +
    xlab("Time Point") +
    ylab("Average log 2Fold Change Expression") +
    ggtitle("Average log 2fc of GO process Genes Over Time") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_x_discrete(labels = c("3 Days", "1 Week", "3 Weeks"))

print(go_genes_plot)




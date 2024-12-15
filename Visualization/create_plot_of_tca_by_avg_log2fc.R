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

# Subset the data to include only the tca genes
tca_genes <- gene_list$tca

# Ensure gene names are uppercase in the dataset
mouse_vcm_all_genes_avg_log2fc$gene <- toupper(mouse_vcm_all_genes_avg_log2fc$gene)

# Filter for tca genes
mouse_vcm_tca_genes_avg_log2fc <- mouse_vcm_all_genes_avg_log2fc %>%
  filter(gene %in% tca_genes)

# Reshape data from wide to long format
mouse_vcm_tca_genes_avg_log2fc_long <- mouse_vcm_tca_genes_avg_log2fc %>%
  pivot_longer(cols = -gene, names_to = "time_point", values_to = "avg_log2fc")

# Set the correct order of time points.
mouse_vcm_tca_genes_avg_log2fc_long$time_point <- factor(mouse_vcm_tca_genes_avg_log2fc_long$time_point,
                                                         levels = c("avg_log2fc_3d",
                                                                    "avg_log2fc_1w",
                                                                    "avg_log2fc_3w"))

# Print to confirm
print(head(mouse_vcm_tca_genes_avg_log2fc_long))

# Create the ggplot for tca genes
ggplot(mouse_vcm_tca_genes_avg_log2fc_long, aes(x = time_point, y = avg_log2fc, group = gene, color = gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set1") +
  xlab("Time Point") +
  ylab("Average log 2Fold Change Expression") +
  ggtitle("Average log 2fc of TCA Genes Over Time") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_discrete(labels = c("3 Days", "1 Week", "3 Weeks"))

# Remember to export the plot
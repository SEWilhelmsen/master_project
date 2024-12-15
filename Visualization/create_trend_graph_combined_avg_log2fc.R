## Create trend graph from combined_avg_log2fc (*ikke ferdig*)
# Input:combined_avg_log2fc data frame
# Output: plot of gene expression over time
# Silje Wilhelmsen

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Import gene list for avg_log2fc
mouse_vcm_all_genes_avg_log2fc <- readRDS("C:/Users/siljeew/RNAseq/mouse_vcm_all_genes_avg_log2fc.rds")


# Specify the gene of interest
gene_of_interest <- "HK1"

# Filter the data for the specific gene
mouse_vcm_one_gene_avg_log2fc <- mouse_vcm_specific_genes_avg_log2fc %>%
  filter(gene == gene_of_interest)

# Print the filtered data to check
print(mouse_vcm_one_gene_avg_log2fc)

# Assuming the column names in your data frame are: "gene", "avg_log2FC_3d", "avg_log2FC_1w", "avg_log2FC_3w"

# Verify column names
colnames(mouse_vcm_one_gene_avg_log2fc)

# Convert data to long format, considering correct column names
mouse_vcm_one_gene_long <- mouse_vcm_one_gene_avg_log2fc %>%
  pivot_longer(cols = starts_with("avg_log2fc_"), names_to = "time_point", values_to = "avg_log2fc")

# Print the long format data to check
print(mouse_vcm_one_gene_long)

# Check that it correctly captured the values and names
print(mouse_vcm_one_gene_long)

# Define the correct time point order
mouse_vcm_one_gene_long$time_point <- factor(mouse_vcm_one_gene_long$time_point, levels = c("avg_log2fc_3d", "avg_log2fc_1w", "avg_log2fc_3w"))

# Arrange the data by time points
mouse_vcm_one_gene_arranged <- mouse_vcm_one_gene_long %>%
  arrange(time_point)

# Print the arranged data to verify
print(mouse_vcm_one_gene_arranged)

# Plot the data using ggplot2
ggplot(mouse_vcm_one_gene_arranged, aes(x = time_point, y = avg_log2fc)) +
  geom_line(group = 1, color = "red", size = 3) +
  geom_point(color = "red", size = 2) +
  labs(title = paste("Trend Graph for Gene:", gene_of_interest),
       x = "Time Point",
       y = "Average Log2 Fold Change") +
  theme_minimal() +
  scale_x_discrete(labels = c("3 Days", "1 Week", "3 Weeks"))
  
## Dont know which of these versions are better
# ggplot(mouse_vcm_one_gene_arranged, aes(x = time_point, y = avg_log2fc, group = gene, color = gene)) +
#   geom_line(size = 2) +
#   geom_point(size = 2) +
#   scale_color_manual(values = c(gene_of_interest = "red")) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#   xlab("Time Point") +
#   ylab("Average log2 Fold Change Expression") +
#   ggtitle("Average log2 Fold Change of Expression over Time") +
#   theme_minimal() +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, size = 1),
#     axis.line = element_line(color = 'black')
#   )

  
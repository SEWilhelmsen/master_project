## Combine significant genes from .xlsx, only the column:  avg_log2FC (*ikke ferdig*)
# Silje Wilhelmsen

library(readxl)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)

## Import data

top_genes_vcm_24h <- read_excel("C:/Users/siljeew/snRNAseq/24h_sham_AB/top_genes_vcm_24h.xlsx")
top_genes_vcm_3d <- read_excel("C:/Users/siljeew/snRNAseq/3days_sham_AB/top_genes_vcm_3d.xlsx")
top_genes_vcm_1w <- read_excel("C:/Users/siljeew/snRNAseq/1week_sham_AB/top_genes_vcm_1w.xlsx")
top_genes_vcm_3w <- read_excel("C:/Users/siljeew/snRNAseq/3weeks_sham_AB/top_genes_vcm_3w.xlsx")


genes_24h <- top_genes_vcm_24h %>%
  dplyr::select(gene, avg_log2FC) 

genes_3d <- top_genes_vcm_3d %>%
  dplyr::select(gene, avg_log2FC)

genes_1w <- top_genes_vcm_1w %>%
  dplyr::select(gene, avg_log2FC)

genes_3w <- top_genes_vcm_3w %>%
  dplyr::select(gene, avg_log2FC)


# Extract gene names into vectors
gene_list_24h <- genes_24h$gene
gene_list_3d <- genes_3d$gene
gene_list_1w <- genes_1w$gene
gene_list_3w <- genes_3w$gene


# Find the intersection of the four gene lists
common_genes <- Reduce(intersect, list(gene_list_24h, gene_list_3d, gene_list_1w,  gene_list_3w))


## Why the FCCCCKKKK funker ikke dette ????!!!!
# Filter the data frames to keep only common genes and rename column to include timepoint
filtered_genes_24h <- genes_24h %>%
  filter(gene %in% common_genes)%>%
  rename(avg_log2fc_24h = avg_log2FC)

filtered_genes_3d <- genes_3d %>%
  filter(gene %in% common_genes)%>%
  rename(avg_log2fc_3d = avg_log2FC)

filtered_genes_1w <- genes_1w %>%
  filter(gene %in% common_genes) %>%
  rename(avg_log2fc_1w = avg_log2FC)

filtered_genes_3w <- genes_3w %>%
  filter(gene %in% common_genes) %>%
  rename(avg_log2fc_3w = avg_log2FC)


#Combine the filtered data frames
merged_genes <- filtered_genes_24h %>%
  inner_join(filtered_genes_3d, by = "gene", suffix = c("_24h", "_3d")) %>%
  inner_join(filtered_genes_1w, by = "gene", suffix = c("", "_1w")) %>%
  inner_join(filtered_genes_3w, by = "gene", suffix = c("", "_3w"))


# Convert the "gene" column to row names
merged_genes <- merged_genes %>%
  column_to_rownames(var = "gene")

# This does not work
#colnames <- c("avg_log2fc_24h", "avg_log2fc_3d", "avg_log2fc_1w", "avg_log2fc_3w")

 
# Save the combined data frame to a CSV file if needed
write.csv(merged_genes, "common_genes_with_avg_log2fc.csv")

#change names of columns to only the time point???

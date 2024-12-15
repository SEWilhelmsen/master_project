## Combine significant genes from .xlsx, only the column:  avg_log2FC
# Silje Wilhelmsen

library(readxl)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)

## Import data

top_genes_VCM_24h <- "C:/Users/siljeew/snRNAseq/24h_sham_AB/top_genes_VCM_24h.xlsx"
top_genes_VCM_3d <- "C:/Users/siljeew/snRNAseq/3days_sham_AB/top_genes_VCM_3d.xlsx"
top_genes_VCM_1w <- "C:/Users/siljeew/snRNAseq/1week_sham_AB/top_genes_VCM_1w.xlsx"
top_genes_VCM_3w <- "C:/Users/siljeew/snRNAseq/3weeks_sham_AB/top_genes_VCM_3w.xlsx"

genes_24h <- read_excel(top_genes_VCM_24h, sheet = 1, col_names = TRUE) %>%
  select(gene, avg_log2FC)

genes_3d <- read_excel(top_genes_VCM_3d, sheet = 1, col_names = TRUE) %>%
  select(gene, avg_log2FC)

genes_1w <- read_excel(top_genes_VCM_1w, sheet = 1, col_names = TRUE) %>%
  select(gene, avg_log2FC)

genes_3w <- read_excel(top_genes_VCM_3w, sheet = 1, col_names = TRUE) %>%
  select(gene, avg_log2FC)

# Extract gene names into vectors
gene_list_24h <- genes_24h$gene
gene_list_3d <- genes_3d$gene
gene_list_1w <- genes_1w$gene
gene_list_3w <- genes_3w$gene

# Find the intersection of the two gene lists
common_genes <- Reduce(intersect, list(gene_list_24h, gene_list_3d, gene_list_1w,  gene_list_3w))


# Filter the data frames to keep only common genes and rename column to include timepoint
filtered_genes_24h <- genes_24h %>%
  filter(gene %in% common_genes)%>%
  rename(avg_log2FC_24h = avg_log2FC)

filtered_genes_3d <- genes_3d %>%
  filter(gene %in% common_genes)%>%
  rename(avg_log2FC_3d = avg_log2FC)

filtered_genes_1w <- genes_1w %>%
  filter(gene %in% common_genes) %>%
  rename(avg_log2FC_1w = avg_log2FC)

filtered_genes_3w <- genes_3w %>%
  filter(gene %in% common_genes) %>%
  rename(avg_log2FC_3w = avg_log2FC)


#Combine the filtered data frames
merged_genes <- filtered_genes_24h %>%
  inner_join(filtered_genes_3d, by = "gene", suffix = c("_24h", "_3d")) %>%
  inner_join(filtered_genes_1w, by = "gene", suffix = c("", "_1w")) %>%
  inner_join(filtered_genes_3w, by = "gene", suffix = c("", "_3w"))

# Convert the "gene" column to row names
merged_genes <- merged_genes %>%
  column_to_rownames(var = "gene")

View(merged_genes)

# Save the combined data frame to a CSV file if needed
write.csv(merged_genes, "common_genes_with_avg_log2FC.csv")

# Print the merged data frame
print(merged_genes)



## Heatmap

# Define GO terms
go_id <- list(
  glycolysis = "GO:0006096",
  FAO = "GO:0019395",
  TCA = "GO:0006099"
)

# Function to retrieve symbols for a GO term
get_genes_for_go <- function(go_id) {
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                      keys = go_id, 
                                      keytype = "GOALL", 
                                      columns = "ENTREZID")[, "ENTREZID"]
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = entrez_ids, 
                                        keytype = "ENTREZID", 
                                        columns = "SYMBOL")[, "SYMBOL"]
  gene_symbols
}

# Retrieve genes for each GO term
go_genes <- lapply(go_id, get_genes_for_go)


# normalize_symbols <- function(genes) {
#   toupper(genes)
# }
# normalized_go_genes <- lapply(go_genes, normalize_symbols)
# 
# all_go_genes <- unique(unlist(normalized_go_genes)) # Unlist GO_genes and create a unique list of GO related genes
# 
# # Function to filter significant genes based on GO terms
# filter_go_genes <- function(merged_genes, all_go_genes) {
#   common_genes <- intersect(normalized_sign_genes, all_go_genes) # Find the intersection between sign_genes and GO_genes
#   sign_genes %>% filter(toupper(gene) %in% common_genes) # Subset sign_genes based on common_genes
# }
# 
# # Ensure row order is consistent
# rownames(filter_go_genes) <- filter_go_genes$gene
# filter_go_genes <- filter_go_genes[ , -1]
# 
# 
# 
# #Find the intersection of the two gene lists
# filter_merged_genes_by_go <- intersect(merged_genes, go_genes)
# 
# 
# # Function to filter merged_genes based on GO term gene lists
# filter_merged_genes_by_go <- function(merged_genes, go_gene_list) {
#   merged_genes %>%
#     filter(rownames(merged_genes) %in% go_gene_list)
# }
# 
# 
# 
# View(filtered_genes) ## Zero rows - where did the rows disappear?
# 
# 
# 
# filter_merged_genes_by_go <- as.matrix(filter_merged_genes_by_go)



row_anno <- rowAnnotation(
  Process = merged_data$process,
  col = list(Process = c("glycolysis" = "coral", "FAO" = "pink", "TCA" = "bisque")),
  show_annotation_name = TRUE
)

heatmap <- Heatmap(filter_merged_genes_by_go,
                   name = "Expression",
                   cluster_rows = FALSE, 
                   show_row_names = TRUE, 
                   row_names_gp = gpar(fontsize = 8),
                   row_names_side = "left",
                   row_title = "Genes",
                   show_row_dend = TRUE,
                   left_annotation = row_anno,
                   cluster_columns = TRUE,  
                   show_column_names = TRUE,
                   column_title = "Expression AB vs SHAM",
                   show_column_dend = FALSE,
                   heatmap_legend_param = list(
                     title = "avg_log2FC",
                     legend_direction = "vertical",
                     legend_height = unit(4, "cm")
                   ),
                   row_gap = unit(2, "mm"),
                   heatmap_height = unit(0.4, "cm")*nrow(filter_merged_genes_by_go),
                   col = colorRamp2(c(min(filter_merged_genes_by_go), 0, max(filter_merged_genes_by_go)),
                                    c("navy", "white", "red"))
)

draw(heatmap)

## Volcano plot
## RNA sequencing data SHAM vs. AB 
## Timepoint: 3 weeks 
## SHAM: 'TP5-1', 'Sample1-TP75B', 'Sample3-TP52A'
## AB: 'Sample2-TP511B', Sample4-TP62A', 'TP51-3'
## Cluster: Ventricular cardiomyocytes
## Silje Wilhelmsen

# Load the data 
load("V:/Silje/snRNAseq/3weeks_sham_AB/Mouse_VCM_3w.h5seurat")

# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(Mouse_VCM_3w) <- Mouse_VCM_3w@meta.data$orig.ident

# Find differentially expressed genes
markers <- FindMarkers(Mouse_VCM_3w, ident.1 = "3 Weeks - AB", ident.2 = "3 Weeks - SHAM")

markers <- markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))

top_genes <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  head(200)


# 0: AB = SHAM, +: mer i AB enn SHAM, -: mindre i AB enn SHAM
ggplot(markers, aes(x = avg_log2FC, y = log_pval)) +
  geom_point(aes(color = p_val_adj < 0.05)) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot VCM - 3 Weeks",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  ylim(0, 50)

png(ggplot, filename = "Vol_SHAM_AB_VCM_3w") #How to export plot?
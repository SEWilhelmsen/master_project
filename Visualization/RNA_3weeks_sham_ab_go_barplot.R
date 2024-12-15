## Gene ontology analysis barplot and dotplot
## RNA sequencing data SHAM vs. AB 
## Timepoint: 3 weeks 
## SHAM: 'TP5-1', 'Sample1-TP75B', 'Sample3-TP52A'
## AB: 'Sample2-TP511B', Sample4-TP62A', 'TP51-3'
## Clustering: Ventricular cardiomyocytes
## Silje Wilhelmsen

# Load the data 
load("C:/Users/siljeew/3weeks_sham_AB/Mouse_VCM_3w.RData")

# Assign the column "orig.ident" as the identity columns for the Seurat object.
Idents(Mouse_VCM_3w) <- Mouse_VCM_3w@meta.data$orig.ident

# Find differentially expressed genes
markers <- FindMarkers(Mouse_VCM_3w, ident.1 = "3 Weeks - AB", ident.2 = "3 Weeks - SHAM")

markers <- markers %>%
  tibble::rownames_to_column(var = "gene") %>%
  mutate(log_pval = -log10(p_val_adj))

markers$RowNames <- rownames(markers)
markers <- FindAllMarkers(Mouse_VCM_3w)
significant_genes <- markers[markers$p_val_adj < 0.05, "gene"]

# GO gene ontology 
gene_list <- bitr(significant_genes, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)

# Perform GO analysis
go_results_VCM <- enrichGO(gene = gene_list$ENTREZID,
                           OrgDb = org.Hs.eg.db,
                           keyType = "ENTREZID",
                           ont = "BP",  # Biological Process
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)

barplot(go_results_VCM,showCategory = 20, title = "Top 20 biological processes in VCM - 3 Weeks")
dotplot(go_results_VCM, showCategory = 20, title = "Top 20 biological processes in VCM - 3 Weeks")

# How to save plots?
##

go_results_df <- as.data.frame(go_results_VCM)
go_results_df_sorted <- go_results_df[order(go_results_df$p.adjust), ]
head(go_results_df_sorted) # This gave the same result as go_results_df. Is this automatically done in the previous code? 

nrow(go_results_df) #Inspect the number of rows/observations in the data frame
ncol(go_results_df) #Inspect the number of columns in the data frame


## Search: "lipid"
lipid_modification_results <- go_results_df[grep("lipid", go_results_df$Description, ignore.case = TRUE), ]

# Sort by p.adjust in descending order
lipid_modification_results <- lipid_modification_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(lipid_modification_results, +
                     aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", 
       y = "p.adjust", 
       title = "Significant Lipid Modification Results within VCM")
# Display the bar plot
print(bar_plot)
# View the filtered results
print(lipid_modification_results)


## Search: "fatty acid"
FA_results <- go_results_df[grep("fatty acid", go_results_df$Description, ignore.case = TRUE), ]
FA_results <- FA_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(FA_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Fatty Acid Results within VCM")
print(bar_plot)
print(FA_results)



## Search: "metabo"
metabo_results <- go_results_df[grep("metabo", go_results_df$Description, ignore.case = TRUE), ]
# Sort the metabolic processes by p.adjust in descending order
metabo_results_sorted <- metabo_results %>%
  mutate(Description = reorder(Description, -p.adjust))
# Create a bar plot with sorted metabolic processes
bar_plot <- ggplot(metabo_results_sorted, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Metabolism Results within VCM")
# Display the bar plot
print(bar_plot)
print(metabo_results)

## Search: "carb"
carb_results <- go_results_df[grep("carb", go_results_df$Description, ignore.case = TRUE), ]
carb_results <- carb_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(carb_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant Carb Results within VCM")
print(bar_plot)
print(carb_results)

## Search: "beta"
beta_results <- go_results_df[grep("beta", go_results_df$Description, ignore.case = TRUE), ]
beta_results <- beta_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(beta_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant beta Results within VCM")
print(bar_plot)
print(beta_results)

## Search: "keto"
keto_results <- go_results_df[grep("keto", go_results_df$Description, ignore.case = TRUE), ]
keto_results <- keto_results %>%
  mutate(Description = reorder(Description, -p.adjust))
bar_plot <- ggplot(keto_results, aes(x = Description, y = p.adjust)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(x = "Metabolism Term", y = "p.adjust", title = "Significant keto Results within VCM")
print(bar_plot)
print(keto_results)

# Create plot of specific genes to combine with plot of go_process
# Silje Wilhelmsen



# Combine the data processing and plotting into a single function for simplicity
create_go_genes_plot <- function(data_avg_logfc, go_process_of_interest_genes) {  
  
  # Subset the data to include only genes in go_process_of_interest_genes
  data_avg_logfc_genes_of_interest <- data_avg_logfc %>%
    filter(gene %in% go_process_of_interest_genes)
  
  # Reshape from wide to long format
  data_avg_logfc_genes_of_interest_go_process_long <- data_avg_logfc_genes_of_interest %>%
    pivot_longer(cols = -gene, names_to = "time_point", values_to = "avg_logfc")
  
  # Set the correct order of timepoints 
  data_avg_logfc_genes_of_interest_go_process_long$time_point <- factor(data_avg_logfc_genes_of_interest_go_process_long$time_point, 
                                                                         levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
  
  
  # Create the ggplot
  go_genes_plot <- ggplot(data_avg_logfc_genes_of_interest_go_process_long, 
                          aes(x = time_point, y = avg_logfc, group = gene, color = gene)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 3) +
    labs(x = NULL,
         y = "Average logFC expression",
         color = "Genes",
         title = "Glycolysis genes") +
    scale_color_brewer(palette = "Set1") +
    scale_color_manual(
        values = c("GCK" = "#E41A1C", "HK1" = "#377EB8", "HK2" = "#4DAF4A", "PFKM" = "#984EA3"),
        labels = expression(italic("Gck"), italic("Hk1"), italic("Hk2"), italic("Pfkm"))) + # Labels in italic
    scale_y_continuous(limits = c(-1, 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 0.5, size = 18, color = "black", margin = margin(t = 15)), 
          axis.text.y = element_text(size = 26, colour = "black"), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 26, margin = margin(r = 15)), # Change margins: r = right, b = bottom, t = top, l = left
          plot.title = element_text(size = 26, margin = margin(b = 20)),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 30),
          axis.line.y.right = element_blank()) 
  
  
  print(go_genes_plot)
}

# Colors:
# scale_color_manual(
#   values = c("CD36" = "#E41A1C"), 
#   labels = expression(italic("Cd36"))) + # Set label to italic
# gene_plot_fatty_acid_uptake <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_fatty_acid_uptake)


# Colors:
# scale_color_manual(
# values = c("CPT1A" = "#E41A1C", "CPT1B" = "#377EB8", "CPT2" = "#4DAF4A"), 
# labels = expression(italic("Cpt1a"), italic("Cpt1b"), italic("Cpt2"))) +  # Labels in italic
# gene_plot_b_oxidation <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_b_oxidation)

# scale_color_manual(
#   values = c("SLC2A4" = "#E41A1C"), 
#   labels = expression(italic("Slc2a4"))) +  # Labels in italic
# gene_plot_glucose_uptake <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_glucose_uptake)


# scale_color_manual(
#   values = c("GCK" = "#E41A1C", "HK1" = "#377EB8", "HK2" = "#4DAF4A", "PFKM" = "#984EA3"), 
#   labels = expression(italic("Gck"), italic("Hk1"), italic("Hk2"), italic("Pfkm"))) + # Labels in italic
gene_plot_glycolysis <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
print(gene_plot_glycolysis)

# scale_color_manual(
#   values = c("DLAT" = "#E41A1C", "DLD" = "#377EB8", "PDHA1" = "#4DAF4A", "PDHX" = "#984EA3"), 
#   labels = expression(italic("Dlat"), italic("Dld"), italic("Pdha1"), italic("Pdhx"))) + # Labels in italic
# gene_plot_pdh <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_pdh)
  

# scale_color_manual(
#   values = c("BDH1" = "#E41A1C"), 
#   labels = expression(italic("Bdh1"))) + # Labels in italic
# gene_plot_bdh1 <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_bdh1)

# scale_color_manual(
#   values = c("OXCT1" = "#E41A1C"), 
#   labels = expression(italic("Oxct1")) + # Labels in italic
# gene_plot_oxct1 <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_oxct1)


# scale_color_manual(
#   values = c("IDH1" = "#E41A1C", "IDH2" = "#377EB8", "IDH3B" = "#4DAF4A"), 
#   labels = expression(italic("Idh1"), italic("Idh2"), italic("Idh3b"))) + # Labels in italic
# gene_plot_tca <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_tca)


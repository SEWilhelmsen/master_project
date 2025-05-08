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
         title = "Tricarboxylic Acid Cycle Genes") +
    scale_color_brewer(palette = "Set1") +
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


# gene_plot_tca <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_tca)

# gene_plot_fatty_acid_uptake <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_fatty_acid_uptake)

# gene_plot_b_oxidation <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_b_oxidation)

# gene_plot_bdh1 <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_bdh1)

# gene_plot_oxct1 <- create_go_genes_plot(data_avg_logfc, go_process_of_interest_genes)
# print(gene_plot_oxct1)

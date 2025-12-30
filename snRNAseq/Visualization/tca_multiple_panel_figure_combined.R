# Create combined plot of tricarboxylic acid cycle 
# Silje Wilhelmsen

# Create barplot of go process grouped by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process.R")
# Create line plot of average 2fold change by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_plot_avglog2fc_for_specific_genes.R")

# Create barplot of go process grouped by stress status:
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process_by_stress.R") #This is not finished
# Create line plot of specific genes grouped by stress status:
## If new gene: 
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_gene_plot_by_stress_.R")
## If updating previous gene:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/tca_multiple_panel_figure.R")



# Make a figure consisting of multiple figures
##############################################################################################
# Set margins
go_plot_tca <- go_plot_tca +
  theme(plot.margin = unit(c(60,20,2,30), "pt"))  # top, right, bottom, left margin
gene_plot_tca <- gene_plot_tca +
  theme(plot.margin = unit(c(60,1,2,110), "pt"))  # top, right, bottom, left margin

single_gene_plot_idh1 <- single_gene_plot_idh1 +
  theme(plot.margin = unit(c(80,1,2,30), "pt"))  # top, right, bottom, left margin
single_gene_plot_idh2 <- single_gene_plot_idh2 +
  theme(plot.margin = unit(c(80,1,2,50), "pt"))  # top, right, bottom, left margin

single_gene_plot_idh3b <- single_gene_plot_idh3b +
  theme(plot.margin = unit(c(80,50,2,30), "pt"))  # top, right, bottom, left margin



# Example of laying out in a 3x3 grid
total_tca_plot <- (go_plot_tca | gene_plot_tca) /
  (single_gene_plot_idh1 | single_gene_plot_idh2) /
  (single_gene_plot_idh3b + plot_spacer())


# Save the combined plot—adjust width and height for space between plots
ggsave(file.path(output_dir_plot, "total_tca_plot_new.pdf"), plot = total_tca_plot, width = 24, height = 24)
ggsave(file.path(output_dir_plot, "total_tca_plot_new.png"), plot = total_tca_plot, width = 24, height = 24)
ggsave(file.path(output_dir_plot, "total_tca_plot_new.tiff"), plot = total_tca_plot, width = 24, height = 24)

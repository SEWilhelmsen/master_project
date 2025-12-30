# Create combined plot of fatty acids
# Silje Wilhelmsen

# Create barplot of go process grouped by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process.R")
# Create line plot of average 2fold change by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_plot_avglogfc_for_specific_genes.R")

# Create barplot of go process grouped by stress status:
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process_by_stress.R") #This is not finished and not used
# Create line plot of specific genes grouped by stress status:
## If new gene: 
# file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_gene_plot_by_stress_.R")
## If updating previous gene:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/pdh_multiple_panel_figure.R")

# Make a figure consisting of multiple figures
##############################################################################################
# Set margins
go_plot_fatty_acid_uptake <- go_plot_fatty_acid_uptake +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
go_plot_b_oxidation <- go_plot_b_oxidation +
  theme(plot.margin = unit(c(1,1,2,5), "cm"))  # top, right, bottom, left margin

gene_plot_fatty_acid_uptake <- gene_plot_fatty_acid_uptake +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
gene_plot_b_oxidation <- gene_plot_b_oxidation +
  theme(plot.margin = unit(c(1,1,2,6), "cm"))  # top, right, bottom, left margin

single_gene_plot_cd36 <- single_gene_plot_cd36 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
single_gene_plot_cpt1a <- single_gene_plot_cpt1a +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin

single_gene_plot_cpt1b <- single_gene_plot_cpt1b +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
single_gene_plot_cpt2 <- single_gene_plot_cpt2 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin


# Example of laying out in a 3x3 grid
total_fatty_acid_plot <- (go_plot_fatty_acid_uptake | go_plot_b_oxidation) /
  (gene_plot_fatty_acid_uptake | gene_plot_b_oxidation) /
  (single_gene_plot_cd36 | single_gene_plot_cpt1a) /
  (single_gene_plot_cpt1b | single_gene_plot_cpt2)



# Save the combined plot—adjust width and height for space between plots
ggsave(file.path(output_dir_plot, "total_fatty_acid_plot_new.png"), plot = total_fatty_acid_plot, width = 24, height = 28)

ggsave(file.path(output_dir_plot, "total_fatty_acid_plot_new.tiff"), plot = total_fatty_acid_plot, width = 24, height = 28)

ggsave(file.path(output_dir_plot, "total_fatty_acid_plot_new.pdf"), plot = total_fatty_acid_plot, width = 24, height = 28)


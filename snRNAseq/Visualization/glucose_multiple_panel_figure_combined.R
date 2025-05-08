# Create combined plot of glucose
# Silje Wilhelmsen

# Create barplot of go process grouped by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_barplot_go_process.R")
# Create line plot of average 2fold change by ORAB vs SHAM:
file.edit("C:/Users/siljeew/Master_project/snRNAseq/Visualization/run_create_plot_avglog2fc_for_specific_genes.R")

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
go_plot_glucose_uptake <- go_plot_glucose_uptake +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
go_plot_glycolysis <- go_plot_glycolysis +
  theme(plot.margin = unit(c(1,1,2,5), "cm"))  # top, right, bottom, left margin

gene_plot_glucose_uptake <- gene_plot_glucose_uptake +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
gene_plot_glycolysis <- gene_plot_glycolysis +
  theme(plot.margin = unit(c(1,1,2,4), "cm"))  # top, right, bottom, left margin

single_gene_plot_slc2a4 <- single_gene_plot_slc2a4 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
single_gene_plot_gck <- single_gene_plot_gck +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin

single_gene_plot_hk1 <- single_gene_plot_hk1 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin
single_gene_plot_hk2 <- single_gene_plot_hk2 +
  theme(plot.margin = unit(c(1,1,2,1), "cm"))  # top, right, bottom, left margin

single_gene_plot_pfkm <- single_gene_plot_pfkm +
  theme(plot.margin = unit(c(1,5,2,1), "cm"))  # top, right, bottom, left margin



# Example of laying out in a 3x3 grid
total_glucose_plot <- (go_plot_glucose_uptake | go_plot_glycolysis) /
  (gene_plot_glucose_uptake | gene_plot_glycolysis) /
  (single_gene_plot_slc2a4 | single_gene_plot_gck) /
  (single_gene_plot_hk1 | single_gene_plot_hk2) /
  (single_gene_plot_pfkm | plot_spacer()) 

# Print the combined plot
print(total_glucose_plot)


# Save the combined plot—adjust width and height for space between plots
ggsave(file.path(output_dir_plot, "total_glucose_plot_new.png"), plot = total_glucose_plot, width = 24, height = 28)

ggsave(file.path(output_dir_plot, "total_glucose_plot_new.pdf"), plot = total_glucose_plot, width = 24, height = 28)
openPDF(file.path(output_dir_plot, "total_glucose_plot_new.pdf"))




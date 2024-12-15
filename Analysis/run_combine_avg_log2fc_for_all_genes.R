# Run combine_avg_log2fc_for_all_genes.R
# Silje Wilhelmsen

# Load libraries
source("load_combine_avg_log2fc_libraries.R")

# Load data and do stuff?
source("combine_avg_log2fc_for_all_genes.R")

# Do analysis?

# Visualize analysis? 
## Plot
source("create_plot_of_glycolysis_by_avg_log2fc.R")
source("create_plot_of_pdh_by_avg_log2fc.R")
source("create_plot_of_tca_by_avg_log2fc.R")
source("create_plot_of_markers_by_avg_log2fc.R")

## Heatmap
source("create_heatmap_combined_avg_log2fc")

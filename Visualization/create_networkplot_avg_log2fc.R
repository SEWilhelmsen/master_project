# Create network plot for genes and process
# Silje Wilhelmsen

# Load libraries 


# Prepare data 

# Create and save plot
##################################
# Function to create plot
create_network <- function() {
  
}

# Create plot
network_plot <- create_network(heatmap_matrix_numeric, go_association_vector)

# Save plot
ggsave(file.path(output_dir_plot, paste(go_processes_of_interest_for_heatmap, "_heatmap.png", sep = "")), plot = combined_plot, width = 10, height = 7)

# Process data for PROGENy analysis
# Silje Wilhelmsen

# Function to process data for each time point
process_timepoint <- function(timepoint_name, layer) {
  # Extract data for current timepoint
  data_matrix <- GetAssayData(
    object = mouse_vcm_all_time_points_with_stress_status, 
    layer = layer,
    assay = "RNA"
  )
  
  # Convert to Dense Matrix
  dense_matrix <- as.matrix(data_matrix)
  
  # Run Progeny Analysis
  pathway_activity_scores <- progeny(dense_matrix, scale = FALSE, organism = "Mouse", top = 500)
  
  # Scale Pathway Activity Scores
  scaled_pathway_scores <- scale(pathway_activity_scores)
  
  # Convert Scores to Data Frame
  progeny_scores_df <- as.data.frame(t(scaled_pathway_scores)) %>%
    tibble::rownames_to_column("Pathway") %>%
    tidyr::gather(Cell, Activity, -Pathway) %>%
    left_join(conditions_df, by = "Cell")  # Join to add condition info to each cell
  
  # Ensure correct grouping and summarization
  summarized_progeny_scores <- progeny_scores_df %>%
    group_by(Pathway, Condition, Timepoint) %>%
    summarise(avg = mean(Activity), std = sd(Activity), .groups = 'drop')
  
  return(progeny_scores_df)
  return(summarized_progeny_scores)
}

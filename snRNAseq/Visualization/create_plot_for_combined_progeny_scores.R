library(dplyr)
library(tidyr)
library(openxlsx)



# # Load data
# mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/Labuser/master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")
# 
# # Define the data layers in a named list
# data_layers <- list(
#   "6 Hours" = "data.mouse_6h_full_merge_project",
#   "12 Hours" = "data.mouse_12h_full_merge_project", 
#   "1 Day" = "data.mouse_12h_full_merge_project",
#   "3 Weeks" = "data.Full Merge Project 3w"
# )
# 
# conditions_df <- data.frame(Cell = colnames(mouse_vcm_all_time_points_with_stress_status),
#                             Stress_Status = mouse_vcm_all_time_points_with_stress_status$Stress_Status,
#                             Timepoint = mouse_vcm_all_time_points_with_stress_status$Timepoint)
# 
# # Function to process data for each time point
# process_timepoint <- function(timepoint_name, layer) {
#   # Extract data for current timepoint
#   data_matrix <- GetAssayData(
#     object = mouse_vcm_all_time_points_with_stress_status, 
#     layer = layer,
#     assay = "RNA"
#   )
#   
#   # Convert to Dense Matrix
#   dense_matrix <- as.matrix(data_matrix)
#   
#   # Run Progeny Analysis
#   pathway_activity_scores <- progeny(dense_matrix, scale = FALSE, organism = "Mouse", top = 500)
#   
#   # Scale Pathway Activity Scores
#   scaled_pathway_scores <- scale(pathway_activity_scores)
#   
#   # Convert Scores to Data Frame
#   progeny_scores_df <- as.data.frame(t(scaled_pathway_scores)) %>%
#     tibble::rownames_to_column("Cell") %>%
#     tidyr::gather(Pathway, Activity, -Cell) %>%
#     inner_join(conditions_df, by = "Cell")  # Join to add condition info to each cell
#   
#   # Group by Pathway and Stress_Status
#   summarized_progeny_scores <- progeny_scores_df %>%
#     group_by(Pathway, Stress_Status, Timepoint) %>%
#     summarise(avg = mean(Activity), std = sd(Activity), .groups = 'drop')
#   
#   # Add Timepoint column based on the current processing timepoint name
#   summarized_progeny_scores_for_xlsx <- summarized_progeny_scores %>%
#     mutate(Timepoint = timepoint_name) %>%
#     select(Timepoint, everything())  # Optionally move Timepoint to the front
#   
#   return(summarized_progeny_scores_for_xlsx)
# }
# 
# # Process each timepoint
# all_results <- lapply(names(data_layers), function(tp) {
#   process_timepoint(tp, data_layers[[tp]])
# }) %>%
#   bind_rows()  # Combine all results into one data frame
# 
# View(all_results)
# 
# # # Write to a single Excel file with all data combined
# # write.xlsx(all_results, "C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_scores_combined.xlsx")













# Prepare data for heatmap
###########################################################
# Remove NA values
non_zero_columns <- combined_progeny_scores %>%
  select(where(~ any(. != 0)))

# Display the modified data frame
head(non_zero_columns)

ggplot(combined_progeny_scores, aes(x = Timepoint, y = Androgen, fill = )) +
  geom_bar(stat = "identity")

# Reshape the data to long format
long_data <- non_zero_columns %>%
  pivot_longer(cols = Androgen:TNFa,  # Specify your range of pathway columns
               names_to = "Pathway",
               values_to = "Value")


# Create heatmap
#######################################################

# Iterate over each unique pathway to create separate heatmaps
unique_pathways <- unique(long_data$Pathway)

# Loop through the pathways
for (pathway in unique_pathways) {
  data_subset <- filter(long_data, Pathway == pathway)
  
  # Create a heatmap for each pathway
  p <- ggplot(data_subset, aes(x = Timepoint, y = Stress_Status, fill = Value)) +
    geom_tile() + # Use geom_tile() for heatmap
    scale_fill_gradient(low = "blue", high = "red") + # Customize colors
    theme_minimal() + # Use a minimal theme
    labs(title = paste("Heatmap for", pathway),
         fill = "Expression") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(p)
}



















# Remove NA values
# Identify and remove columns that sum to zero across all rows

# Assuming combined_data includes the necessary columns
combined_data$Timepoint <- factor(combined_data$Timepoint, 
                                  levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))


# Verify the levels are set correctly
levels(combined_data$Timepoint)

# Remove NA values
# Identify and remove columns that sum to zero across all rows
non_zero_cols <- colSums(combined_data != 0) > 0
filtered_df <- combined_data[, non_zero_cols]

# Print the resulting data frame
print(filtered_df)

# Check if all columns contain numeric data
str(filtered_df)


# Create heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(
  seq(min(combined_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(combined_data, na.rm = TRUE) / paletteLength, max(combined_data, na.rm = TRUE), 
      length.out = floor(paletteLength / 2))
)


split = mtcars$cyl

pheatmap(
  combined_data, 
  cluster_rows = FALSE,  # Preserve order of Stress_Status
  cluster_cols = FALSE,  # The order here is defined by the levels in Timepoint
  split = combined_data,
  fontsize = 14, 
  fontsize_row = 10, 
  color = myColor,
  breaks = progenyBreaks, 
  main = "Androgen Pathway Activity by Stress Status and Time", 
  angle_col = 45
)







# Pathway:Androgen
# Make sure 'Set' column is a factor in the correct time order
library(dplyr)
library(tidyr)
library(pheatmap)

# Assuming combined_data includes the necessary columns
combined_data$Timepoint <- factor(combined_data$Timepoint, 
                                  levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Verify the levels are set correctly
levels(combined_data$Timepoint)

androgen_data <- combined_data %>%
  select(Stress_Status, Timepoint, Androgen) %>%
  pivot_wider(names_from = Timepoint, values_from = Androgen) %>%
  column_to_rownames(var = "Stress_Status")


paletteLength <- 100
myColor <- colorRampPalette(c("darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(
  seq(min(androgen_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
  seq(max(androgen_data, na.rm = TRUE) / paletteLength, max(androgen_data, na.rm = TRUE), 
      length.out = floor(paletteLength / 2))
)

# Create heatmap
p1 <- pheatmap(
  androgen_data, 
  cluster_rows = FALSE,  # Preserve order of Stress_Status
  cluster_cols = FALSE,  # The order here is defined by the levels in Timepoint
  fontsize = 14, 
  fontsize_row = 10, 
  color = myColor,
  breaks = progenyBreaks, 
  main = "Androgen Pathway Activity by Stress Status and Time", 
  angle_col = 45
)

# Render the plot
print(p1)





# Pathway: EGFR
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Reshape the data so Stress_Status is rows and Set is columns
egfr_data <- combined_data %>%
  select(Stress_Status, Timepoint, EGFR) %>%
  pivot_wider(names_from = Timepoint, values_from = EGFR) %>%
  column_to_rownames(var = "Stress_Status")

# Check structure to verify numeric columns
str(egfr_data)

# Define color palette for heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(egfr_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
                   seq(max(egfr_data, na.rm = TRUE) / paletteLength, 
                       max(egfr_data, na.rm = TRUE), 
                       length.out = floor(paletteLength / 2)))

# Create heatmap
p2 <- pheatmap(egfr_data, 
         cluster_rows = FALSE,  # Keep original order of Stress_Status
         cluster_cols = FALSE,  # Keep original order of Set
         fontsize = 14, 
         fontsize_row = 10, 
         color = myColor,
         breaks = progenyBreaks, 
         main = "EGFR Pathway Activity by Stress Status and Time", 
         angle_col = 45)

p2





# Pathway: JAK-STAT
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Reshape the data so Stress_Status is rows and Set is columns
jakstat_data <- combined_data %>%
  select(Stress_Status, Timepoint, 'JAK-STAT') %>%
  pivot_wider(names_from = Timepoint, values_from = 'JAK-STAT') %>%
  column_to_rownames(var = "Stress_Status")

# Check structure to verify numeric columns
str(jakstat_data)

# Define color palette for heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(jakstat_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
                   seq(max(jakstat_data, na.rm = TRUE) / paletteLength, 
                       max(jakstat_data, na.rm = TRUE), 
                       length.out = floor(paletteLength / 2)))

# Create heatmap
p3 <- pheatmap(jakstat_data, 
               cluster_rows = FALSE,  # Keep original order of Stress_Status
               cluster_cols = FALSE,  # Keep original order of Set
               fontsize = 14, 
               fontsize_row = 10, 
               color = myColor,
               breaks = progenyBreaks, 
               main = "JAK-STAT Pathway Activity by Stress Status and Time", 
               angle_col = 45)

p3





# Pathway: NFkB
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))


# Reshape the data so Stress_Status is rows and Set is columns
nfkb_data <- combined_data %>%
  select(Stress_Status, Timepoint, NFkB) %>%
  pivot_wider(names_from = Timepoint, values_from = NFkB) %>%
  column_to_rownames(var = "Stress_Status")

# Check structure to verify numeric columns
str(nfkb_data)

# Define color palette for heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(nfkb_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
                   seq(max(nfkb_data, na.rm = TRUE) / paletteLength, 
                       max(nfkb_data, na.rm = TRUE), 
                       length.out = floor(paletteLength / 2)))

# Create heatmap
p4 <- pheatmap(nfkb_data, 
               cluster_rows = FALSE,  # Keep original order of Stress_Status
               cluster_cols = FALSE,  # Keep original order of Set
               fontsize = 14, 
               fontsize_row = 10, 
               color = myColor,
               breaks = progenyBreaks, 
               main = "NFkB Pathway Activity by Stress Status and Time", 
               angle_col = 45)

p4




# Pathway: TNFa
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))


# Reshape the data so Stress_Status is rows and Set is columns
tnfa_data <- combined_data %>%
  select(Stress_Status, Timepoint, TNFa) %>%
  pivot_wider(names_from = Timepoint, values_from = TNFa) %>%
  column_to_rownames(var = "Stress_Status")

# Check structure to verify numeric columns
str(tnfa_data)

# Define color palette for heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(tnfa_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
                   seq(max(tnfa_data, na.rm = TRUE) / paletteLength, 
                       max(tnfa_data, na.rm = TRUE), 
                       length.out = floor(paletteLength / 2)))

# Create heatmap
p5 <- pheatmap(tnfa_data, 
               cluster_rows = FALSE,  # Keep original order of Stress_Status
               cluster_cols = FALSE,  # Keep original order of Set
               fontsize = 14, 
               fontsize_row = 10, 
               color = myColor,
               breaks = progenyBreaks, 
               main = "TNFa Pathway Activity by Stress Status and Time", 
               angle_col = 45)

p5






# Pathway: p53
combined_data$Timepoint <- factor(combined_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

# Reshape the data so Stress_Status is rows and Set is columns
p53_data <- combined_data %>%
  select(Stress_Status, Timepoint, p53) %>%
  pivot_wider(names_from = Timepoint, values_from = p53) %>%
  column_to_rownames(var = "Stress_Status")

# Check structure to verify numeric columns
str(p53_data)

# Define color palette for heatmap
paletteLength <- 100
myColor <- colorRampPalette(c("Darkblue", "white", "red"))(paletteLength)
progenyBreaks <- c(seq(min(p53_data, na.rm = TRUE), 0, length.out = ceiling(paletteLength / 2) + 1),
                   seq(max(p53_data, na.rm = TRUE) / paletteLength, 
                       max(p53_data, na.rm = TRUE), 
                       length.out = floor(paletteLength / 2)))

# Create heatmap
p6 <- pheatmap(p53_data, 
               cluster_rows = FALSE,  # Keep original order of Stress_Status
               cluster_cols = FALSE,  # Keep original order of Set
               fontsize = 14, 
               fontsize_row = 10, 
               color = myColor,
               breaks = progenyBreaks, 
               main = "p53 Pathway Activity by Stress Status and Time", 
               angle_col = 45)

p6



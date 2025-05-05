# ANOVA + TukeyHSD test for progeny scores for one specific pathway
# Silje Wilhelmsen

# Load libraries
library(dplyr)
library(writexl)
library(readxl)



# Perform ANOVA + TukeyHSD
#########################################################################################
# Perform ANOVA+TukeyHSD test without prior summarization
perform_anova_tukey_timepoint_focus <- function(data, group_col, activity_col, time_col, pathways_of_interest) {
  anova_results <- list()
  tukey_results <- list()
  
  for (pathway in pathways_of_interest) {
    # Subset data for the specific pathway of interest
    pathway_data <- data[data$Pathway == pathway, ]
    pathway_data[[group_col]] <- factor(pathway_data[[group_col]])  # Ensure Stress_Status is a factor
    
    # Split data by the given timepoint
    split_data <- split(pathway_data, pathway_data[[time_col]])
    path_tukey_results <- list()
    
    for (time_point in names(split_data)) {
      data_subset <- split_data[[time_point]]
      
      print(paste("Processing Timepoint:", time_point))
      
      if (nrow(data_subset) >= 3) {  # Ensure enough data for ANOVA
        # Perform ANOVA, focusing on Stress_Status
        formula <- as.formula(paste(activity_col, "~", group_col))
        anova_result <- aov(formula, data = data_subset)
        print(summary(anova_result))
        
        # Tukey HSD post-hoc testing on Stress_Status
        tukey_result <- TukeyHSD(anova_result)
        
        if ("Stress_Status" %in% names(tukey_result)) {
          tukey_df <- as.data.frame(tukey_result[["Stress_Status"]])
          tukey_df$comparison <- rownames(tukey_df)
          
          # Add group1 and group2 from comparison parts
          comparison_parts <- strsplit(tukey_df$comparison, "-")
          tukey_df$group1 <- sapply(comparison_parts, `[`, 1)
          tukey_df$group2 <- sapply(comparison_parts, `[`, 2)
          
          # Add pathway and timepoint to results
          tukey_df$Pathway <- pathway
          tukey_df$Timepoint <- time_point
          
          # Collect results
          path_tukey_results[[time_point]] <- tukey_df
        }
      }
    }
    
    # Bind results by timepoint for each pathway and collect
    path_tukey_df <- bind_rows(path_tukey_results)
    tukey_results[[length(tukey_results) + 1]] <- path_tukey_df
  }
  
  tukey_results_df <- bind_rows(tukey_results)
  
  # Save the results to file
  write_xlsx(tukey_results_df, tukey_xlsx_path)
  
  return(tukey_results_df)
}


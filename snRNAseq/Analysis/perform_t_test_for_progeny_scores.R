# t-test for progeny scores for one specifict pathway

# Load libraries
library(dplyr)
library(writexl)
library(readxl)

# Define parameters
#########################################################################################
# Define pathways
pathways_of_interest <- c("p53", "TNFa", "NFkB", "JAK-STAT", "EGFR", "Androgen")

# Define groups to compare
sham_group <- "SHAM - CM"
stressed_group <- "Stressed CM"


# Perform t-test
#########################################################################################
test_results <- list()

# Perform t-test on the two groups
for (pathway in pathways_of_interest) {
  data_for_t_test <- progeny_scores_df %>%
    filter(Stress_Status %in% c(sham_group, stressed_group)) %>%
    filter(Pathway == pathway) %>%
    group_by(Pathway, Timepoint) %>%
    filter(n_distinct(Stress_Status) == 2) %>%
    ungroup()
  
  test_result <- data_for_t_test %>%
    group_by(Timepoint) %>%
    summarise(
      t_test = list(t.test(Activity ~ Stress_Status, data = .)),
      mean_sham = mean(Activity[Stress_Status == sham_group]),
      mean_stressed = mean(Activity[Stress_Status == stressed_group]),
      .groups = 'drop'
    ) 
  
  for (i in seq_along(test_result$t_test)) {
    p_value <- test_result$t_test[[i]]$p.value
    mean_sham <- test_result$mean_sham[i]
    mean_stressed <- test_result$mean_stressed[i]
    timepoint <- test_result$Timepoint[i]
    
    test_results[[length(test_results) + 1]] <- list(
      Pathway = pathway,
      Timepoint = timepoint,
      P_Value = p_value,
      Mean_Sham = mean_sham,
      Mean_Stressed = mean_stressed,
      Groups_Compared = paste(sham_group, "vs", stressed_group)
    )
  }
}

test_results_df <- do.call(rbind, lapply(test_results, as.data.frame))
test_results_df$P_Adjusted <- p.adjust(test_results_df$P_Value, method = "BH")

# Save results
#############################################################################
file_path <- "C:/Users/siljeew/Master_project/snRNAseq/Data/progeny_scores_t_test.xlsx"

if (file.exists(file_path)) {
  # Read existing results
  existing_results_df <- read_excel(file_path)
  
  # Ensure column data types match
  existing_results_df <- existing_results_df %>%
    mutate(
      Pathway = as.character(Pathway),
      Timepoint = as.character(Timepoint)
    )
  
  # Update results
  updated_results_df <- existing_results_df %>%
    full_join(test_results_df, by = c("Pathway", "Timepoint", "Groups_Compared"), suffix = c("", ".new")) %>%
    mutate(
      P_Value = ifelse(is.na(P_Value.new), P_Value, P_Value.new),
      P_Adjusted = ifelse(is.na(P_Adjusted.new), P_Adjusted, P_Adjusted.new),
      Mean_Sham = ifelse(is.na(Mean_Sham.new), Mean_Sham, Mean_Sham.new),
      Mean_Stressed = ifelse(is.na(Mean_Stressed.new), Mean_Stressed, Mean_Stressed.new)
    ) %>%
    select(-ends_with(".new"))  # Remove the '.new' columns
} else {
  updated_results_df <- test_results_df  # No existing file, use the test results
}

# Save the updated results
write_xlsx(updated_results_df, file_path)

cat("Results have been saved and/or updated in progeny_scores_t_test.xlsx\n")
<<<<<<< HEAD
# t-test for progeny scores for one specific pathway
# Silje Wilhelmsen

# group1 <- "SHAM CM"
# group2 <- "Not stressed CM"
# group2 <- "Stressed CM"

# Prepare a list to store the results
t_test_results <- list()

# Perform t-tests for each pathway at each timepoint
for (pathway in unique(progeny_results_filtered$Pathway)) {
  # Subset data for the current pathway
  pathway_data <- progeny_results_filtered %>%
    filter(Pathway == pathway, Stress_Status %in% c(group1, group2))
  
  # Split by Timepoint and perform t-tests
  for (timepoint in unique(pathway_data$Timepoint)) {
    data_subset <- pathway_data %>%
      filter(Timepoint == timepoint)
    
    # Ensure there are enough observations per group
    group1_values <- data_subset %>% filter(Stress_Status == group1) %>% pull(Activity)
    group2_values <- data_subset %>% filter(Stress_Status == group2) %>% pull(Activity)
    
    if (length(group1_values) >= 2 && length(group2_values) >= 2) {
      # Perform the t-test
      t_test_result <- t.test(group2_values, group1_values)
      
      # Store results
      t_test_results[[length(t_test_results) + 1]] <- data.frame(
        Pathway = pathway,
        Timepoint = timepoint,
        p_value = t_test_result$p.value,
        t_statistic = t_test_result$statistic,
        group2 = group2,
        mean_group2 = mean(group2_values, na.rm = TRUE),
        group1 = group1,
        mean_group1 = mean(group1_values, na.rm = TRUE),
        conf_low = t_test_result$conf.int[1],
        conf_high = t_test_result$conf.int[2],
        stringsAsFactors = FALSE
      )
    } else {
      message(paste("Not enough data for pathway:", pathway, "at timepoint:", timepoint))
    }
  }
}

# Combine all results into a single data frame
t_test_results_df <- bind_rows(t_test_results)

# Perform Benjamini-Hochberg correction for p-values
if (nrow(t_test_results_df) > 0) {
  t_test_results_df$p_adjusted_bh <- p.adjust(t_test_results_df$p_value, method = "BH")
}

View(t_test_results_df)
# View results
print(t_test_results_df)


# Add label
t_test_results_df <- t_test_results_df %>%
  mutate(p_adjusted_label_bh = ifelse(p_adjusted_bh < 0.001, "***",
                                      ifelse(p_adjusted_bh < 0.01, "**",
                                             ifelse(p_adjusted_bh < 0.05, "*",
                                                    "ns"))))


=======
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
>>>>>>> 7df848ad59037d1ccdbeb933fda380e63ccfea21

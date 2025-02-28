# Combine PROGENy scores. 
# The PROGENy scores are weighted scores based on predefined weigths of genes. Each cell gets a score that is then summarized. 
# The scores are compared to the means of the whole dataset. 
# PROGENy scores above 0 are increased compared to the mean of the entire dataset, and below 0 are lower than the mean. 

library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# # Combine the data frames by row
# combined_data <- rbind(progeny_scores_12h, progeny_scores_6h, progeny_scores_1d, progeny_scores_3d, progeny_scores_1w, progeny_scores_3w)
# 
# # Save the combined data to an Excel sheet
# # Create a new workbook
# workbook <- createWorkbook()
# 
# # Add a sheet and write the combined data
# addWorksheet(workbook, "Combined_Progeny_Scores")
# writeData(workbook, sheet = "Combined_Progeny_Scores", combined_data)
# 
# # Save the workbook to an xlsx file
# saveWorkbook(workbook, file = "C:/Users/Labuser/snRNAseq/Data/combined_progeny_scores.xlsx", overwrite = TRUE)
# 
# view(combined_data)

# Read the combined data. 
combined_progeny_scores <- read_excel("C:/Users/siljeew/Master_project/snRNAseq/Data/combined_progeny_scores.xlsx")
View(combined_progeny_scores)

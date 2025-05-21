# Create plot of distribution of stress status in SHAM
# Silje Wilhelmsen

# Load libraries
library(ggplot2)
library(dplyr)
library(Seurat)


# Load data
mouse_vcm_all_time_points_with_stress_status <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_time_points_with_stress_status.Rds")
SeuratObj.all <- mouse_vcm_all_time_points_with_stress_status 


# Add stress status
####################################################################################
unique(SeuratObj.all@meta.data$Group)

# Define genes to indicate stress and thresholds
genes_of_interest <- c("MYH7", "NPPA", "NPPB", "ANKRD1")
VentricularCM <- subset(SeuratObj.all, subset = predicted.celltype.l2 == "Ventricular Cardiomycoyte")
thresholds <- apply(FetchData(VentricularCM, vars = genes_of_interest), 2, function(x) mean(x, na.rm = TRUE))
selected_nuclei <- WhichCells(VentricularCM, expression = MYH7 > 5* thresholds["MYH7"] |
                                NPPA > 5*  thresholds["NPPA"] |
                                NPPB >5*  thresholds["NPPB"] |
                                ANKRD1 > 5* thresholds["ANKRD1"])
meta_data <- SeuratObj.all@meta.data

# Ensure the necessary columns exist
if (!all(c("predicted.celltype.l2", "Group") %in% colnames(meta_data))) {
  stop("Required columns not found in metadata")
}

# Create a subset of SHAM nuclei
sham_subset <- subset(SeuratObj.all, subset = Group == "SHAM")

# Assign 'Stressed CM' if criteria are met
meta_data$Stress_Status[meta_data$predicted.celltype.l2 == "Ventricular Cardiomycoyte" & 
                          rownames(meta_data) %in% selected_nuclei & 
                          meta_data$Group == "SHAM"] <- "Stressed CM"

# Assign 'Non-CM' if the nucleus is not a Ventricular Cardiomyocyte
meta_data$Stress_Status[meta_data$predicted.celltype.l2 != "Ventricular Cardiomycoyte"] <- "Non-CM"

# Assign 'Not stressed CM' if it's a Ventricular Cardiomyocyte but not in the selected nuclei (Group AB)
meta_data$Stress_Status[meta_data$predicted.celltype.l2 == "Ventricular Cardiomycoyte" & 
                          !(rownames(meta_data) %in% selected_nuclei) & 
                          meta_data$Group == "SHAM"] <- "Not stressed CM"


# Assign updated metadata back to the Seurat object
sham_subset@meta.data <- meta_data

# VentricularCM <- subset(sham_subset, subset = predicted.celltype.l2 == "Ventricular Cardiomycoyte")

# stressed_cm_seurat <- subset(sham_subset, cells = selected_nuclei)




# Create plot
#####################################################################################
sham_percentage_data <- sham_subset@meta.data %>%
  filter(Group == "SHAM") %>%
  group_by(Timepoint, Stress_Status) %>%
  summarise(Cell_Count = n()) %>%
  group_by(Timepoint) %>%
  mutate(Total_Cells = sum(Cell_Count),
         Percentage_of_total_cells = (Cell_Count / Total_Cells) * 100)

print(sham_percentage_data)

# sham_percentage_data <- sham_percentage_data %>%
#   mutate(Stress_Status = ifelse(Stress_Status == "SHAM - CM", "SHAM CM", Stress_Status))
# View(sham_percentage_data)

# Set order of stress status
sham_percentage_data$Stress_Status <- factor(sham_percentage_data$Stress_Status,
                                        levels = c("Not stressed CM", "Stressed CM"))

sham_percentage_data$Timepoint <- factor(sham_percentage_data$Timepoint, levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))

stress_in_sham_plot <- ggplot(sham_percentage_data, aes(x = Timepoint, y = Percentage_of_total_cells, fill = Stress_Status)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Distribution of stress status in SHAM",
       x = NULL,
       y = "Stress status (%)",
       fill = "Stress status") +
  scale_y_continuous(breaks = seq(0, 100, by = 10), labels = paste0(seq(0, 100, by = 10), "%")) +
  scale_fill_manual(values = c("Not stressed CM" = "slategray2", "Stressed CM" = "slategray4")) +
  theme(axis.text.x = element_text(hjust = 0.5, size = 26, color = "black", angle = 30, vjust = 0.8, margin = margin(t = 15)),
        axis.text.y = element_text(size = 26, colour = "black"),
        axis.title.y = element_text(size = 26, colour = "black", margin = margin(r = 15)),
        plot.title = element_text(size = 30, color = "black", margin = margin(b = 20, t = 10)),
        legend.title = element_text(size = 30, colour = "black"),
        legend.text = element_text(size = 30, color = "black"),
        axis.line.y.right = element_blank())

print(stress_in_sham_plot)

ggsave(file.path(output_dir_plot, paste("distribution_of_stress_in_sham.png", sep = "")), plot = stress_in_sham_plot, width = 12, height = 10, dpi = 400)
ggsave(file.path(output_dir_plot, paste("distribution_of_stress_in_sham.tiff", sep = "")), plot = stress_in_sham_plot, width = 12, height = 10, dpi = 400)


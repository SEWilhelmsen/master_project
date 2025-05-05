# Create histograms for cell size data
# Input: dataset of small size (e.g. n=10 per group)
# Output: 
# Silje Wilhelmsen


# Preparation
###############################################################################
# file.edit("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/prepare_cell_size_data.R")

cell_size_data_filtered <- read.xlsx("C:/Users/Labuser/master_project/Phenotypic_data/cell_size_data_filtered.xlsx")
head(cell_size_data_filtered)




# Area
##########################################################################################
cell_size_data_filtered$`area_um2` <- as.numeric(cell_size_data_filtered$`area_um2`)

timepoint_order <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")

cell_size_data_for_plot <- cell_size_data_filtered %>%
  filter(area_um2 < 20000) %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoint_order))


# Create histogram of area
histogram <- ggplot(cell_size_data_for_plot, aes(x = area_um2, fill = Condition)) +
  geom_histogram(binwidth = 500, color = "black", alpha = 0.7, position = "identity") +
  labs(title = "Distribution of Cell Areas by Condition at 1 Week",
       x = "Area (um2)",
       y = "Frequency") +
  facet_wrap(~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c( "SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

histogram

ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_size_distribution.pdf")
ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_size_distribution.png")




# Length
##########################################################################################
cell_size_data_filtered$`length_um` <- as.numeric(cell_size_data_filtered$`length_um`)

timepoint_order <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")

cell_size_data_for_plot <- cell_size_data_filtered %>%
  filter(length_um < 1000) %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoint_order))

# Create histogram of area
histogram <- ggplot(cell_size_data_for_plot, aes(x = length_um, fill = Condition)) +
  geom_histogram(binwidth = 20, color = "black", alpha = 0.7, position = "identity") +
  labs(title = "Distribution of Cell Length",
       x = "Length (um2)",
       y = "Frequency") +
  facet_wrap(~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c( "SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

histogram

ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_length_distribution.pdf")
ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_length_distribution.png")




# Width
##########################################################################################
cell_size_data_filtered$`width_um` <- as.numeric(cell_size_data_filtered$`width_um`)

timepoint_order <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")

cell_size_data_for_plot <- cell_size_data_filtered %>%
  filter(width_um < 200) %>%
  mutate(Timepoint = factor(Timepoint, levels = timepoint_order))


# Create histogram of area
histogram <- ggplot(cell_size_data_for_plot, aes(x = width_um, fill = Condition)) +
  geom_histogram(binwidth = 5, color = "black", alpha = 0.7, position = "identity") +
  labs(title = "Distribution of Cell Width",
       x = "Width (um2)",
       y = "Frequency") +
  facet_wrap(~ Timepoint, ncol = 6) +
  scale_fill_manual(values = c( "SHAM" = "grey22", "ORAB" = "coral")) +
  theme_minimal()

histogram

ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_width_distribution.pdf", dpi = 400)
ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/cell_width_distribution.png", dpi = 400)

### QQ-plot to investigate distribution of a variable
# Input: 
# Output: qq-plot of variable
## Silje Wilhelmsen

# Load libraries
library(dplyr)
library(ggpubr)
library(ggplot2)


filtered_data <- processed_data %>%
  filter(condition == "sham" & time_point == "3_weeks") %>%
  select(all_of(c("sample_id", "condition", "time_point", "hw")))
View(filtered_data)

## qq-plot
qqplot_hw <- ggqqplot(filtered_data, 
                      x = "hw", 
                      ggtheme = theme_bw()) + 
  ggtitle("Sham at 3_weeks, HW (mg)") +
  labs(x = "Theoretical Quantiles", 
       y = "Sample Quantiles")
print(qqplot_hw)

output_file <- "Y:/Silje/snRNAseq/Plots/Distribution/qqplot_sham_3_weeks_hw.png"

# Save
ggsave(filename = output_file, plot = qqplot_hw, device = "png", width = 20, height = 20, units = "cm")


## histogram

histogram_hw <- ggplot(filtered_data, aes(x = hw)) +
  geom_histogram(binwidth = 5, colour = "black", fill = "grey") + 
  ggtitle("Sham at 3_weeks, HW (mg)") +
  labs(x = "hw",
       y = "Frequency") +
  theme_bw()

print(histogram_hw)

# Define the file path where the histogram will be saved
output_file <- "Y:/Silje/snRNAseq/Plots/Distribution/histogram_sham_3_weeks_hw.png"

# Save the histogram
ggsave(filename = output_file, plot = histogram_hw, device = "png", width = 20, height = 20, units = "cm")





######################################################
# Welch test for unequal variance between the groups #
######################################################


# Filter datasets for condition "sham" and "AB"
data_sham <- processed_data %>% 
  filter(condition == "sham" & time_point == "1_week") %>%
  select(sample_id, condition, time_point, hw)

data_ab <- processed_data %>% 
  filter(condition == "AB" & time_point == "1_week") %>%
  select(sample_id, condition, time_point, hw)

# Perform the t-test
t.test(data_sham$hw, data_ab$hw, var.equal = FALSE)


# Filter datasets for condition "sham" and "AB"
data_sham <- processed_data %>% 
  filter(condition == "sham" & time_point == "3_weeks") %>%
  select(sample_id, condition, time_point, hw)

data_ab <- processed_data %>% 
  filter(condition == "AB" & time_point == "3_weeks") %>%
  select(sample_id, condition, time_point, hw)


# Perform the t-test
t.test(data_sham$hw, data_ab$hw, var.equal = FALSE)



############################################################
# T-test
###########################################################



# Filter datasets for condition "sham" and "AB"
data_sham <- processed_data %>% 
  filter(condition == "sham" & time_point == "6_hours") %>%
  select(sample_id, condition, time_point, hw)

data_ab <- processed_data %>% 
  filter(condition == "AB" & time_point == "6_hours") %>%
  select(sample_id, condition, time_point, hw)

# Perform the t-test
t.test(data_sham$hw, data_ab$hw, var.equal = TRUE)


# Filter datasets for condition "sham" and "AB"
data_sham <- processed_data %>% 
  filter(condition == "sham" & time_point == "12_hours") %>%
  select(sample_id, condition, time_point, hw)

data_ab <- processed_data %>% 
  filter(condition == "AB" & time_point == "12_hours") %>%
  select(sample_id, condition, time_point, hw)


# Perform the t-test
t.test(data_sham$hw, data_ab$hw, var.equal = TRUE)



#########################################################
#                   Boxplot 
########################################################


# Do test and store into something

p_values <- data.frame(
  time_point = c("6_hours", "12_hours", "1_day", "3_days", "1_week", "3_weeks"),
  group1 = "sham", 
  group2 = "AB",
  p = c(0.1542, 0.8254, 0.001681, 0.002226, 0.0000448, 0.00001046)
)

p_values <- p_values %>%
  mutate(y.position = c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6))

View(p_values)

# # Arrange the results to get the correct positions for the significances
# stat.test <- stat.test %>% 
#   add_xy_position(x = "condition", dodge = 0.8)


highlighted_samples <- processed_data %>%
  filter(Used_already_for == "Multiome seq") %>%
  pull(sample_id)

processed_data_highlighted <- processed_data %>%
  mutate(is_highlighted = ifelse(sample_id %in% highlighted_samples, TRUE, FALSE))

View(processed_data_highlighted)

# Create the boxplot with significance annotations
boxplot_hw <- ggplot(processed_data_highlighted, aes(x = condition, y = hw, fill = condition)) +
  geom_boxplot(alpha = 0.7, show.legend = FALSE) +
  geom_point(data = filter(processed_data_highlighted, is_highlighted),
             aes(x = condition, y = hw),
             size = 2, shape = 21, fill = "red") +
  facet_grid(. ~ time_point) +
  labs(title = "Heart weight (mg) ratio of all samples", 
       x = "condition", 
       y = "hw") +
  scale_fill_manual(name = "condition", values = c("sham" = "white", "AB" = "black"), 
                    labels = c("sham" = "Sham", "AB" = "AB")) +
  theme_minimal() +
  theme(legend.position = "right") +
  stat_pvalue_manual(p_values, label = "p", tip.length = 0.01)

# Display the plot
print(boxplot_hw)



# Basic boxplot
simple_boxplot <- ggplot(processed_data_highlighted, aes(x = condition, y = hw, fill = condition)) +
  geom_boxplot(alpha = 0.7) +
  facet_grid(.~time_point)
  

# Display the simple plot
print(simple_boxplot)


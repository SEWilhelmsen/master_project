# Create boxplots for cell size data
# Load prepared cell size data, 
# Silje Wilhelmsen

# Preparation
#######################################################
# Load libraries
source("C:/Users/siljeew/Master_project/Phenotypic_data/Analysis/load_libraries_for_plots_phenotype_data.R")
load_libraries()

# Read data
cell_size_data_filtered <- read_excel("C:/Users/siljeew/Master_project/Phenotypic_data/Data/cell_size_data_filtered.xlsx")
head(cell_size_data_filtered)

# Set factors and orders
cell_size_data_filtered$Timepoint <- factor(cell_size_data_filtered$Timepoint, 
                                            levels = c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks"))
cell_size_data_filtered$Condition <- factor(cell_size_data_filtered$Condition, 
                                            levels = c("SHAM", "ORAB"))


# Create boxplot for area
#########################################################################################
# Make sure data is numeric
cell_size_data_filtered$`area_um2` <- as.numeric(cell_size_data_filtered$`area_um2`)

# Manually set significance label
p_adj_label <- data.frame(
  Timepoint = c("3 Days", "3 Weeks"),
  label = c("***", "*"),
  Condition = c("SHAM")
)

# Jitter plot
jitter_plot_area <- ggplot(cell_size_data_filtered, aes(x = Timepoint, y = area_um2, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 1, show.legend = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.9), 
              size = 0.6, alpha = 0.6,
              show.legend = FALSE) +
  labs(title = "Cell Area", 
       x = NULL, 
       y = "Area (\U00B5m\U00B2)") +
  scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
                    labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
  theme_classic() +
  theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
        strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
        panel.border = element_blank(),
        axis.title.y = element_text(size = 26, color = "black", margin = margin(t = 15)),
        axis.text.x = element_text(hjust = 1, size = 26, color = "black", angle = 30, margin = margin(t = 5)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 26, colour = "black", margin = margin(r = 5)),
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 30), 
        axis.line.y.right = element_blank()) +
  geom_text(data = p_adj_label, aes(y = 10000, x = Timepoint, label = label), 
            size = 16, vjust = 0.8, color = "black")


# Display the plot
print(jitter_plot_area)

ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/jitter_cell_area.png", plot = jitter_plot_area)


# boxplot <- ggplot(cell_size_data_for_plot, aes(x = Timepoint, y = area_um2, fill = Condition)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   #geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
#               #color = "black", size = 0.4, alpha = 0.6) +
#   labs(title = "Cell size", 
#        x = "Timepoint", 
#        y = "Area (µm {supsc('2')})") +
#   scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
#                     labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#   theme_classic() +
#   theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
#         strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
#         panel.border = element_blank(),
#         axis.title.x = element_text(size = 26, color = "black", margin = margin(t = 15)),
#         axis.text.x = element_text(hjust = 0.5, size = 20, color = "black", margin = margin(t = 5)),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(size = 28, colour = "black", margin = margin(r = 5)),
#         axis.line.y.right = element_blank())
# 
# # Display the plot
# print(boxplot)
# 
# #Save plot 
# ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/box_cell_area.png", width = 19, height = 12, dpi = 400)
# ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/box_cell_area.pdf", width = 19, height = 12, dpi = 400)
# 
# # Violin plot
# violin_plot <- ggplot(cell_size_data_filtered, aes(x = Timepoint, y = area_um2, fill = Condition)) +
#   geom_violin(linewidth = 1) + 
#   geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.5) +
#   
#     labs(title = "Cell Area", 
#        x = NULL, 
#        y = "Area (µm\U00B2)") +
#   scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
#                     labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#   theme_classic() +
#   theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
#         strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
#         panel.border = element_blank(),
#         axis.title.y = element_text(size = 26, color = "black", margin = margin(t = 15)),
#         axis.text.x = element_text(hjust = 0.5, size = 20, color = "black", margin = margin(t = 5)),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 28, colour = "black", margin = margin(r = 5)),
#         legend.title = element_text(size = 30), 
#         legend.text = element_text(size = 30), 
#         axis.line.y.right = element_blank()) +
#   geom_text(data = p_adj_label, aes(y = 10000, x = Timepoint, label = label), 
#             size = 18, vjust = 0.8, color = "black") 
# 
# 
# # Display the plot
# print(violin_plot)
# 
# ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/violin_cell_area.png", plot = violin_plot)






# Width
################################################################################################
cell_size_data_filtered$width_um <- as.numeric(cell_size_data_filtered$width_um)

# Manually set significance label
p_adj_label <- data.frame(
  Timepoint = c("12 Hours", "3 Days"),
  label = c("***", "**"),
  Condition = c("SHAM")
)

# Jitter plot
jitter_plot_width <- ggplot(cell_size_data_filtered, aes(x = Timepoint, y = width_um, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 1) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.9), 
              size = 0.6, alpha = 0.6,
              show.legend = FALSE) +
  labs(title = "Cell Width", 
       x = NULL, 
       y = "Width (\U00B5m)") +
  scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
                    labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
  theme_classic() +
  theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
        strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
        panel.border = element_blank(),
        axis.title.y = element_text(size = 26, color = "black", margin = margin(t = 15)),
        axis.text.x = element_text(hjust = 1, size = 26, color = "black", angle = 30, margin = margin(t = 5)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 26, colour = "black", margin = margin(r = 5)),
        legend.title = element_text(size = 40), 
        legend.text = element_text(size = 40),
        legend.position = "bottom",
        axis.line.y.right = element_blank()) +
  geom_text(data = p_adj_label, aes(y = 100, x = Timepoint, label = label), 
            size = 16, vjust = 0.8, color = "black")

# Display the plot
print(jitter_plot_width)


# boxplot <- ggplot(cell_size_data_for_plot, aes(x = Timepoint, y = width_um, fill = Condition)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
#               color = "black", size = 0.4, alpha = 0.6) +
#   
#   labs(title = "Cell size", 
#        x = "Timepoint", 
#        y = "width_um") +
#   
#   scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
#                     labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#   
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"),
#     plot.title = element_text(size = 16, face = "bold")
#   ) +
#   stat_compare_means(aes(group = Condition), method = "t.test", label = "p.format")
# 
# # Display the plot
# print(boxplot)
# 
# #Save plot 
# ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/box_cell_width.png", width = 19, height = 12)
# ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/box_cell_width.pdf", width = 19, height = 12)





# Length
#############################################################################################
cell_size_data_filtered$`length_um` <- as.numeric(cell_size_data_filtered$`length_um`)

# Manually set significance label
p_adj_label <- data.frame(
  Timepoint = c("12 Hours", "1 Day", "3 Days", "3 Weeks"),
  label = c("*", "**", "***", "***"),
  Condition = c("SHAM")
)

jitter_plot_length <- ggplot(cell_size_data_filtered, aes(x = Timepoint, y = length_um, fill = Condition)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 1, show.legend = FALSE) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.7, dodge.width = 0.9), 
              size = 0.6, alpha = 0.6,
              show.legend = FALSE) +
  labs(title = "Cell Length", 
       x = NULL, 
       y = "Length (\U00B5m)") +
  scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
                    labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
  theme_classic() +
  theme(plot.title = element_text(size = 30, margin = margin(b = 15, t = 10), hjust = 0.5),
        strip.text = element_text(size = 30, color = "black"),  # Adjust timepoint title
        panel.border = element_blank(),
        axis.title.y = element_text(size = 26, color = "black", margin = margin(t = 15)),
        axis.text.x = element_text(hjust = 1, size = 26, color = "black", angle = 30, margin = margin(t = 5)),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 26, colour = "black", margin = margin(r = 5)),
        legend.title = element_text(size = 30), 
        legend.text = element_text(size = 30), 
        axis.line.y.right = element_blank()) +
  geom_text(data = p_adj_label, aes(y = 350, x = Timepoint, label = label), 
            size = 16, vjust = 0.8, color = "black")


# Display the plot
print(jitter_plot_length)


# boxplot <- ggplot(cell_size_data_for_plot, aes(x = Timepoint, y = length_um, fill = Condition)) +
#   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
#   geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), 
#               color = "black", size = 0.4, alpha = 0.6) +
#   
#   labs(title = "Cell size", 
#        x = "Timepoint", 
#        y = "length_um") +
#   
#   scale_fill_manual(name = "Condition", values = c("SHAM" = "grey22", "ORAB" = "coral"),
#                     labels = c("SHAM" = "SHAM", "ORAB" = "ORAB")) +
#   
#   theme_minimal() +
#   theme(
#     legend.position = "right",
#     axis.text.x = element_text(size = 12),
#     axis.text.y = element_text(size = 12),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"),
#     plot.title = element_text(size = 16, face = "bold")
#   ) +
#   stat_compare_means(aes(group = Condition), method = "t.test", label = "p.format")
# 
# # Display the plot
# print(boxplot)
# 
# #Save plot 
# ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/box_cell_length.png", width = 19, height = 12)
# ggsave("C:/Users/Labuser/master_project/Phenotypic_data/Plots/box_cell_length.pdf", width = 19, height = 12)
# 



# Combine jitter for area, width and length
##################################################################################
jitter_plot_area <- jitter_plot_area +
  theme(plot.margin = unit(c(1,1,3,1), "cm"))  # top, right, bottom, left margin

jitter_plot_width <- jitter_plot_width +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) # This plot has smaller margin because of the legends

jitter_plot_length <- jitter_plot_length +
  theme(plot.margin = unit(c(1,1,3,1), "cm")) 


combined_plot <- plot_grid(jitter_plot_area, 
                           jitter_plot_width, 
                           jitter_plot_length, 
                           ncol = 3)
# combined_plot

ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/cell_size_combined_jitter.png", plot = combined_plot, width = 26, height = 14, dpi = 300)
ggsave("C:/Users/siljeew/Master_project/Phenotypic_data/Plots/cell_size_combined_jitter.tiff", plot = combined_plot, width = 26, height = 14, dpi = 300)


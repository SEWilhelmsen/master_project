# Create plot of estimated proportion of a celltype in sectors

# Important to run processing of data first: 
# open: "C:/Users/Labuser/master_project/Spatials/run_preprocessing_of_spatial_slides.R". From here you get the seurat object "object"
# Preprocessing: "C:/Users/Labuser/master_project/Spatials/run_preprocessing_of_spatial_slides.R")
# Make sure the plot has the right rotation: C:/Users/Labuser/master_project/Spatials/Analysis/rotation_of_slide.R
# Define sectors: C:/Users/Labuser/master_project/Spatials/Analysis/define_sectors.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)


# Load data
# RCTD <- readRDS("C:/Users/Labuser/master_project/Spatials/tmp/A1_RCTD_temp.rds")
load("C:/Users/Labuser/master_project/SCreference_l2_ALLNUCLEI_10112025.RData")
# RCTD <- TP39.4_RCTD_temp

# Make sure you have these objects before continuing:
barcodes <- colnames(RCTD@spatialRNA@counts) # list of spatial barcodes
weights <- RCTD@results$weights # Weights for each cell type per barcode
norm_weights <- normalize_weights(weights)
cell_type_names <- colnames(norm_weights) # List of cell types
subset_df <- as.data.frame(t(as.data.frame(norm_weights[1:2,])))
subset_df$celltypes <- rownames(subset_df); rownames(subset_df) <- NULL



# Verify that norm_weights has 1 as sum
range(norm_weights, na.rm = TRUE)
row_sums_norm <- rowSums(norm_weights, na.rm = TRUE)
summary(row_sums_norm)
table(round(row_sums_norm, 3))
print(colnames(norm_weights)) # Print available celltypes

# Define the celltype of interest
celltype <- "Capillary Endothelial"

if (!celltype %in% colnames(norm_weights)) {
  stop("Could not find celltype in norm_weights: ", celltype, "Available celltypes are: ", paste(colnames(norm_weights), collapse = ", "))
}

# Create a data frame of spots and the celltype proportion
df_prop <- data.frame(spot = rownames(norm_weights),
                      prop = as.numeric(norm_weights[, celltype]),
                      stringsAsFactors = FALSE)

meta <- object@meta.data
meta$spot <- rownames(meta)

df <- df_prop %>%
  left_join(meta %>% select(spot, SectorZone), by = "spot")

n_unmatched <- sum(is.na(df$SectorZone))
if (n_unmatched > 0) {
  warning(n_unmatched, " spots mangler SectorZone (kan skyldes mismatch i spot-navn).")
}

  # Summarize per SectorZone
# prop = estimert andel av spot-signalet som modellen tilskriver den celletype
summary_df <- df %>%
  group_by(SectorZone) %>%
  summarise(
    n_spots = n(),
    n_valid = sum(!is.na(prop)),
    mean_prop = ifelse(n_valid>0, mean(prop, na.rm = TRUE), NA_real_),
    median_prop = ifelse(n_valid>0, median(prop, na.rm = TRUE), NA_real_),
    sd_prop = ifelse(n_valid>0, sd(prop, na.rm = TRUE), NA_real_),
    pct_majority = ifelse(n_valid>0, 100 * sum(prop > 0.5, na.rm = TRUE) / n_valid, NA_real_),
    pct_present10 = ifelse(n_valid>0, 100 * sum(prop > 0.10, na.rm = TRUE) / n_valid, NA_real_)
  ) %>%
  arrange(desc(mean_prop))
print(summary_df)


# Create a data frame of one point per SectorZone 
label_df <- coords %>%
  group_by(SectorZone) %>%
  summarise(x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            Sector = first(Sector),
            Zone = first(Zone),
            n_spots = n()) %>%
  ungroup()


label_df <- label_df %>%
  left_join(summary_df %>% select(SectorZone, mean_prop), by = "SectorZone")


# Format labels
label_df <- label_df %>%
  mutate(label = ifelse(is.na(mean_prop), NA_character_,
                        paste0(round(mean_prop * 100, 1), "%")))

# Create plot
pB_labels <- ggplot(coords, aes(x = x, y = y, color = SectorZone)) +
  geom_point(size = 1.2, alpha = 0.9) +
  geom_path(data = circle_df, aes(x = x, y = y, group = radius),
            color = "black", size = 0.5, linetype = "dashed") +
  annotate("point", x = xc, y = yc, colour = "red", size = 2) +
  coord_equal() +
  facet_wrap(~ Zone) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle(paste0("Estimated proportion of ", celltype))

pB_labels

# Add labels to plot
pB_labels <- pB_labels +
  geom_label(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4,
    color = "black",
    fill = "white",
    alpha = 0.8,
    label.r = unit(0.12, "lines"),
    label.padding = unit(0.15, "lines")
  )

pB_labels + scale_color_viridis_d(option = "turbo")   # eller "magma", "viridis"



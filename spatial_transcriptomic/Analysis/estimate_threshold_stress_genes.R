# Save threshold values as .rds

data_dir <- "C:/Users/siljeew/Master_project/spatial_transcriptomic/Data/"


# bygge en data.frame med ønsket kolonne-rekkefølge
object_df <- data.frame(
  sample              = unique(object@meta.data$sample),
  group               = unique(object@meta.data$group),
  timepoint           = unique(object@meta.data$Timepoint),
  threshold_nppa      = threshold_nppa,
  threshold_nppb      = threshold_nppb,
  threshold_ankrd1    = threshold_ankrd1,
  threshold_myh7      = threshold_myh7,
  stringsAsFactors = FALSE
)

print(object_df)

thresholds_stress <- readRDS("C:/Users/siljeew/Master_project/spatial_transcriptomic/Data/thresholds_stress.rds")

# thresholds_df <- data.frame(
#   sample            = character(),
#   group             = character(),
#   timepoint         = character(),
#   threshold_nppa    = numeric(),
#   threshold_nppb    = numeric(),
#   threshold_ankrd1  = numeric(),
#   threshold_myh7    = numeric(),
#   stringsAsFactors = FALSE
# )

thresholds_stress <- rbind(thresholds_stress, object_df)

print(thresholds_stress)

saveRDS(thresholds_stress, file = paste0(data_dir, "thresholds_stress.rds"))


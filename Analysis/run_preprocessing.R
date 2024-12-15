## Run preprocessing (main)

source("load_libraries_preprocessing.R")
source("preprocessing.R")
source("3_weeks.R")

# "Config"
options(future.globals.maxSize = 5 * 1024^3)  # Increase memory limit to 5 GB

# Create and save Seurat objects in the list
seurat_object_per_mouse <- create_mouse_seurat_object(experiment_data, time_point)

# # Verify the Seurat objects are in the list
# print(seurat_object_per_mouse$TP51)

# Annotation of nuclei - Assignment to each nuclei a score on what type of cell type the nucleus comes from
predictions_list <- annotation_of_nuclei_mouse(seurat_object_per_mouse)
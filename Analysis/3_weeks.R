## List of list
## Timepoint: 3 weeks

## SHAM: 'TP5-1', 'Sample1-TP75B', 'Sample3-TP52A'
## AB: 'Sample2-TP511B', Sample4-TP62A', 'TP51-3'
# Silje Wilhelmsen

# Create a list of lists
experiment_data <- list(
  TP51 = list(id = "TP51", path = "W:/CardioTarget/Multiome seq data/TP5-1/outs/filtered_feature_bc_matrix.h5", condition = "SHAM")
  # TP75 = list(id = "TP75", path = "W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample1-TP75B/outs/filtered_feature_bc_matrix.h5", condition = "SHAM"),
  # TP52 = list(id = "TP52", path = "W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample3-TP52A/outs/filtered_feature_bc_matrix.h5", condition = "SHAM"),
  # TP511 = list(id = "TP511", path = "W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample2-TP511B/outs/filtered_feature_bc_matrix.h5", condition = "AB"),
  # TP62 = list(id = "TP62", path = "W:/CardioTarget/Multiome seq data/GCF0742_GCF0795_CellRanger/Sample4-TP62A/outs/filtered_feature_bc_matrix.h5", condition = "AB"),
  # TP513 = list(id = "TP513", path = "W:/CardioTarget/Multiome seq data/TP51-3/outs/filtered_feature_bc_matrix.h5", condition = "AB")
)

time_point <- "3 weeks"
# Perform UMAP based on subset of seurat data

# Make sure you have defined sectors and zones here:
# C:/Users/Labuser/master_project/Spatials/Analysis/define_sectors.R

print(unique(object@meta.data$SectorZone))
# [1] NA                            "Septum_Outer"                "Septum_Mid"                  "Mid anterior_Mid"           
# [5] "Mid anterior septum_Mid"     "Mid inferior_Mid"            "Mid post_Mid"                "Mid post_Central"           
# [9] "Mid inferior_Outer"          "Septum_Central"              "Mid lateral_Mid"             "Mid anterior septum_Outer"  
# [13] "Mid anterior_Central"        "Mid lateral_Central"         "Mid inferior_Central"        "Mid anterior septum_Central"
# [17] "Mid anterior_Outer"  

seurat_obj <- subset(object, subset = SectorZone == "Mid anterior_Mid")
seurat_obj <- subset(object, subset = SectorZone == "Mid anterior septum_Mid")
seurat_obj <- subset(object, subset = SectorZone == "Septum_Mid")

DefaultAssay(seurat_obj) <- "Spatial.008um"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, reduction.name = "pca.008um") # Time: 4+ hr 
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca.008um", dims = 1:20, graph.name = "SNN.008um", k.param = 30) # Time: 20 min. dims = 1:30 produced error.

names(seurat_obj@graphs)

seurat_obj <- FindClusters(seurat_obj, graph.name = "SNN.008um", resolution = 0.6, cluster.name = "seurat_cluster.008um") # Time: 8 hrs
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca.008um", reduction.name = "umap.008um", dims = 1:20) # Time: 4 hrs


dim.plot <- DimPlot(seurat_obj, reduction = "umap.008um", group.by = "seurat_cluster.008um", label = TRUE,
                    repel = T) + NoLegend()
dim.plot

cluster.plot <- SpatialDimPlot(seurat_obj, group.by = "seurat_cluster.008um", label = FALSE, pt.size.factor = 2) +
  theme(legend.position = "right")

cluster.plot




# Make p5 and p6 for 2-3 clusters from the violin plot, and save individual images if possible instead of displaying all in one plot
# can highlight by identity for more detailed visualization
p0 <- SpatialDimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, expression = seurat_cluster.008um ==
                                                            0), pt.size.factor = 2) + NoLegend() + ggtitle("Cluster 0")
p0

p1 <- SpatialDimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, expression = seurat_cluster.008um ==
                                                            1), pt.size.factor = 2) + NoLegend() + ggtitle("Cluster 1")
p1
p2 <- SpatialDimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, expression = seurat_cluster.008um ==
                                                            2), pt.size.factor = 2) + NoLegend() + ggtitle("Cluster 2")
p2
p3 <- SpatialDimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, expression = seurat_cluster.008um ==
                                                            3), pt.size.factor = 2) + NoLegend() + ggtitle("Cluster 3")
p3
p4 <- SpatialDimPlot(seurat_obj, cells.highlight = WhichCells(seurat_obj, expression = seurat_cluster.008um ==
                                                            8), pt.size.factor = 2) + NoLegend() + ggtitle("Cluster 8")
p4


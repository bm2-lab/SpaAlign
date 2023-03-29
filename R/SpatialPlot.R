SpatialPlot <-
function (initial_clustering_result, sample_information_coordinate) {
    sce_seurat <- initial_clustering_result$sce_seurat
    Idents(sce_seurat) <- initial_clustering_result$sample_information
    require(ggplot2)
    sce_seurat <- RunUMAP(sce_seurat, dims = 1:15)
    sample_information_coordinate$cluster <- Idents(sce_seurat)
    cluster_plot <- DimPlot(sce_seurat, reduction = "umap")
    coordinate_plot <- ggplot(sample_information_coordinate, 
        aes(x = X, y = Y, colour = cluster)) + geom_point(size = 1.2)
    plot(cluster_plot)
    plot(coordinate_plot)
}

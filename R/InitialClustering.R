InitialClustering <-
function (expression_profile, user_offered = FALSE, sample_information_user_offered=NULL, nfeatures = 2000, 
    resolution = (ifelse(user_offered==FALSE,2,FALSE))) 
{
    initial_clustering_result <- list()
    initial_clustering_result$expression_profile <- expression_profile
    require(Seurat)
    sce_seurat <- CreateSeuratObject(expression_profile)
    sce_seurat <- FindVariableFeatures(sce_seurat, nfeatures = nfeatures)
    sce_seurat <- ScaleData(sce_seurat)
    sce_seurat <- RunPCA(sce_seurat, features = VariableFeatures(object = sce_seurat), 
        ndims.print = 1, nfeatures.print = 1)
    if (user_offered) {
        initial_clustering_result$sample_information <- sample_information_user_offered
    }
    else {
        sce_seurat <- FindNeighbors(sce_seurat, reduction = "pca")
        sce_seurat <- FindClusters(sce_seurat, resolution = resolution)
        sample_cluster <- Idents(sce_seurat)
        names(sample_cluster) <- colnames(expression_profile)
        initial_clustering_result$sample_information <- sample_cluster
    }
    initial_clustering_result$high_var_genes <- VariableFeatures(object = sce_seurat)
    initial_clustering_result$sce_seurat <- sce_seurat
    return(initial_clustering_result)
}

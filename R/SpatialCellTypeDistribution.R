SpatialCellTypeDistribution <-
function (sample_information_coordinate, sequence_resolution = c("single_cell","spot"), 
    sample_information_cellType = NULL, sample_information_decon = NULL, 
    neighbour_search_method = (ifelse(sequence_resolution == 
        "single_cell", "KNN", "radius")), k = (ifelse(neighbour_search_method == 
        "KNN", 30, FALSE)), r = (ifelse(neighbour_search_method == 
        "radius", 2, FALSE))) 
{
    SpatialKNN <- function(sample_information_coordinate) {
        require(Seurat)
        test_coordinate <- sample_information_coordinate
        test_coordinate_expand <- test_coordinate
        for (i in 1:5) {
            test_coordinate_expand <- cbind(test_coordinate, 
                test_coordinate_expand)
        }
        colnames(test_coordinate_expand) <- paste(colnames(test_coordinate_expand), 
            1:ncol(test_coordinate_expand))
        sce_seurat <- CreateSeuratObject(t(test_coordinate_expand))
        sce_seurat <- ScaleData(sce_seurat)
        sce_seurat <- RunPCA(sce_seurat, features = rownames(sce_seurat), 
            ndims.print = 1, nfeatures.print = 1, npcs = 5)
        sce_seurat@reductions$pca@cell.embeddings <- as.matrix(test_coordinate)
        sce_seurat <- FindNeighbors(sce_seurat, reduction = "pca", 
            dims = 1:2, return.neighbor = T, k.param = 50)
        knn_sample <- sce_seurat@neighbors$RNA.nn@nn.idx
        knn_value <- sce_seurat@neighbors$RNA.nn@nn.dist
        cell_names <- sce_seurat@neighbors$RNA.nn@cell.names
        row.names(knn_sample) <- cell_names
        row.names(knn_value) <- cell_names
        knn_sample <- t(apply(knn_sample, 1, function(x, y) {
            y <- cell_names
            return(y[x])
        }))
        return(list(knn_value = knn_value, knn_sample = knn_sample))
    }
    SpatialKNN_result <- SpatialKNN(sample_information_coordinate)
    require(reshape2)
    neighbour_search_method = neighbour_search_method[1]
    knn_sample <- SpatialKNN_result$knn_sample
    knn_value <- SpatialKNN_result$knn_value
    knn_choose <- list()
    if (neighbour_search_method == "KNN") {
        if (k < 2) {
            print("k must be bigger than 1")
            break
        }
        for (i in 1:nrow(knn_sample)) {
            knn_choose[[i]] <- knn_sample[i, 1:k]
        }
    }
    else if (neighbour_search_method == "radius") {
        for (i in 1:nrow(knn_sample)) {
            knn_choose[[i]] <- knn_sample[i, knn_value[i, ] <= 
                r]
        }
    }
    else {
        print("Parameter 'neighbour_search_method' must be 'KNN' or 'radius'!")
        break
    }
    names(knn_choose) <- row.names(knn_sample)
    knn_choose_cellType <- knn_choose
    sequence_resolution <- sequence_resolution[1]
    if (sequence_resolution == "single_cell") {
        names(knn_choose_cellType) <- sample_information_cellType[names(knn_choose)]
        for (i in 1:length(knn_choose_cellType)) {
            knn_choose_cellType[[i]] <- sample_information_cellType[knn_choose[[i]]]
        }
        cellType_names <- names(table(sample_information_cellType))
        knn_matrix_list <- list()
        knn_matrix_cellType_list <- list()
        for (i in 1:length(cellType_names)) {
            knn_matrix_cellType <- knn_choose_cellType[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix <- knn_choose[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix_list[[i]] <- knn_matrix
            knn_matrix_cellType_list[[i]] <- knn_matrix_cellType
        }
        names(knn_matrix_list) <- cellType_names
        names(knn_matrix_cellType_list) <- cellType_names
        knn_matrix <- list(knn_matrix_cellType_list = knn_matrix_cellType_list, knn_matrix_list = knn_matrix_list)

        cellType_names <- names(table(sample_information_cellType))
        cell_type_distribution_each_list <- list()
        cell_type_distribution <- matrix(NA, 0, length(cellType_names))
        colnames(cell_type_distribution) <- cellType_names
        cell_type_distribution_rownames <- c()
        for (k in 1:length(cellType_names)) {
            test_cellType_matrix <- knn_matrix[[1]][[k]]
            test_matrix <- knn_matrix[[2]][[k]]
            cell_type_distribution_rownames <- c(cell_type_distribution_rownames, 
                names(test_matrix))
            cell_type_distribution_each <- matrix(0, length(test_cellType_matrix), 
                length(cellType_names))
            colnames(cell_type_distribution_each) <- cellType_names
            for (i in 1:length(test_cellType_matrix)) {
                cellType_table <- table(factor(test_cellType_matrix[[i]], 
                  levels = cellType_names))
                if (sum(cellType_table) == 0) {
                  cell_type_distribution_each[i, ] <- cellType_table
                }
                else {
                  cell_type_distribution_each[i, ] <- cellType_table/sum(cellType_table)
                }
            }
            row.names(cell_type_distribution_each) <- names(test_matrix)
            cell_type_distribution_each_list[[k]] <- cell_type_distribution_each
            cell_type_distribution <- rbind(cell_type_distribution, 
                cell_type_distribution_each)
        }
        names(cell_type_distribution_each_list) <- cellType_names
        row.names(cell_type_distribution) <- cell_type_distribution_rownames
        return(cell_type_distribution)
    }
    else if(sequence_resolution=="spot"){
        sample_information_cellType<-rep(c("a","b","c"),length.out = nrow(sample_information_decon))
        names(sample_information_cellType)<-row.names(sample_information_decon)
        names(knn_choose_cellType) <- sample_information_cellType[names(knn_choose)]
        for (i in 1:length(knn_choose_cellType)) {
            knn_choose_cellType[[i]] <- sample_information_cellType[knn_choose[[i]]]
        }
        cellType_names <- names(table(sample_information_cellType))
        knn_matrix_list <- list()
        knn_matrix_cellType_list <- list()
        for (i in 1:length(cellType_names)) {
            knn_matrix_cellType <- knn_choose_cellType[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix <- knn_choose[names(knn_choose_cellType) == 
                cellType_names[i]]
            knn_matrix_list[[i]] <- knn_matrix
            knn_matrix_cellType_list[[i]] <- knn_matrix_cellType
        }
        names(knn_matrix_list) <- cellType_names
        names(knn_matrix_cellType_list) <- cellType_names
        knn_matrix <- list(knn_matrix_cellType_list = knn_matrix_cellType_list, 
            knn_matrix_list = knn_matrix_list)

        cellType_names <- colnames(sample_information_decon)
        cell_type_distribution <- matrix(NA, 0, length(cellType_names))
        colnames(cell_type_distribution) <- cellType_names
        cell_type_distribution_rownames <- c()
        for (k in 1:length(knn_matrix[[2]])) {
            test_matrix <- knn_matrix[[2]][[k]]
            cell_type_distribution_rownames <- c(cell_type_distribution_rownames, 
                names(test_matrix))
            for (i in 1:length(test_matrix)) {
                if (length(test_matrix[[i]]) == 1) {
                  cell_type_distribution <- rbind(cell_type_distribution, 
                    sample_information_decon[test_matrix[[i]], 
                      ])
                }
                else {
                  cell_type_distribution <- rbind(cell_type_distribution, 
                    apply(sample_information_decon[test_matrix[[i]], 
                      ], 2, mean))
                }
            }
        }
        row.names(cell_type_distribution) <- cell_type_distribution_rownames
        return(cell_type_distribution)
    }
    else{
        print("sequence_resolution must be one of 'single_cell' and 'spot'!")
        break
    }
}

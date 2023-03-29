SpatialCellTypeDistribution_multiple <-
function (sample_information_coordinate_list, sequence_resolution = c("single_cell","spot"), 
    sample_information_cellType_list = NULL, sample_information_decon_list = NULL) 
{
    sequence_resolution<-sequence_resolution[1]
    cell_type_distribution_list<-list()
    if(sequence_resolution=='single_cell'){
        if(is.null(sample_information_cellType_list)){
            print("'sample_information_cellType_list' shouldn't be NULL!")
        }else{
            for(i in 1:length(sample_information_coordinate_list)){
                cell_type_distribution_list[[i]]<-SpatialCellTypeDistribution(sample_information_coordinate = sample_information_coordinate_list[[i]],
                                                                              sequence_resolution = 'single_cell',sample_information_cellType  = sample_information_cellType_list[[i]]) 
            }
        }
    }else if(sequence_resolution=='spot'){
        if(is.null(sample_information_decon_list)){
            print("'sample_information_decon_list' shouldn't be NULL!")
        }else{
            for(i in 1:length(sample_information_coordinate_list)){
                cell_type_distribution_list[[i]]<-SpatialCellTypeDistribution(sample_information_coordinate = sample_information_coordinate_list[[i]],
                                                                              sequence_resolution = 'spot',sample_information_decon  = sample_information_decon_list[[i]]) 
            }
        }
        
    }
    names(cell_type_distribution_list)<-names(sample_information_coordinate_list)
    cellType_union <- colnames(cell_type_distribution_list[[1]])
    for (i in 2:length(cell_type_distribution_list)) {
        cellType_union <- union(cellType_union, colnames(cell_type_distribution_list[[i]]))
    }
    for (j in 1:length(cell_type_distribution_list)) {
        if (length(cellType_union) - length(colnames(cell_type_distribution_list[[j]])) > 
            0) {
            missing_data_1 <- matrix(0, nrow(cell_type_distribution_list[[j]]), 
                length(cellType_union) - length(colnames(cell_type_distribution_list[[j]])))
            colnames(missing_data_1) <- setdiff(cellType_union, 
                colnames(cell_type_distribution_list[[j]]))
            cell_type_distribution_list[[j]] <- cbind(cell_type_distribution_list[[j]], 
                missing_data_1)
            cell_type_distribution_list[[j]] <- cell_type_distribution_list[[j]][, 
                cellType_union]
        }
        row.names(cell_type_distribution_list[[j]])<-paste(names(cell_type_distribution_list)[j],row.names(cell_type_distribution_list[[j]]),sep="_")
    }
    cell_type_distribution_combine <- cell_type_distribution_list[[1]]
    datasets_name <- names(cell_type_distribution_list)
    datasets_lable <- c(rep(datasets_name[1], nrow(cell_type_distribution_list[[1]])))
    names(datasets_lable) <- row.names(cell_type_distribution_list[[1]])
    for (k in 2:length(cell_type_distribution_list)) {
        cell_type_distribution_combine <- rbind(cell_type_distribution_combine, 
            cell_type_distribution_list[[k]])
        datasets_lable2 <- c(rep(datasets_name[k], nrow(cell_type_distribution_list[[k]])))
        names(datasets_lable2) <- row.names(cell_type_distribution_list[[k]])
        datasets_lable <- c(datasets_lable, datasets_lable2)
    }
    return(list("cell_type_distribution_combine" = cell_type_distribution_combine, 
        "datasets_lable"= factor(datasets_lable)))
}

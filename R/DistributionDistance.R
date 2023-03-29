DistributionDistance <-
function (cell_type_distribution, distance = c("jensen-shannon", "manhattan")) 
{
    method_choose <- distance[1]
    ### high accuracy, but lower speed
    if (method_choose == "jensen-shannon") {
        require(philentropy)
        propor_dis <- philentropy::distance(cell_type_distribution, method = method_choose)
        row.names(propor_dis) <- row.names(cell_type_distribution)
        colnames(propor_dis) <- row.names(cell_type_distribution)
        return(propor_dis)
    }
    ### high speed, but lower accuracy
    else if (method_choose == "manhattan") {
        propor_dis <- dist(x = cell_type_distribution, method = "manhattan")
        propor_dis <- as.matrix(propor_dis)
        return(propor_dis)
    }
    else {
        print("method must be jensen-shannon or manhattan")
        break
    }
}

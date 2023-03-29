DomainHclust <-
function (distribution_distance, autoselection = TRUE, auto_resolution=c("low","median","high"),domain_num = ifelse(autoselection = FALSE, 
    10, FALSE)) 
{
    result_hc <- hclust(d = as.dist(distribution_distance), method = "ward.D2")
    if (autoselection) {
        require(dynamicTreeCut)
        auto_resolution<-auto_resolution[1]
        if(auto_resolution=="low"){
            cluster_hc_dy <- cutreeDynamic(result_hc, distM = distribution_distance,deepSplit = 0)   
        }else if(auto_resolution=="median"){
            cluster_hc_dy <- cutreeDynamic(result_hc, distM = distribution_distance,deepSplit = 1)   
        }else if(auto_resolution=="high"){
            cluster_hc_dy <- cutreeDynamic(result_hc, distM = distribution_distance,deepSplit = 2)   
        }else{
            break
            print("auto_resolution must one of 'low','median','high'.")
        }
        cluster_hc_dy_num <- length(table(cluster_hc_dy))
    }
    else {
        cluster_hc_dy_num = domain_num
    }
    cluster_hc_all = cutree(result_hc, k = 2:cluster_hc_dy_num)
    cluster_hc_all_order <- cluster_hc_all[result_hc$order, ]
    cluster_hc_all_order_rename <- as.data.frame(cluster_hc_all_order)
    colnames(cluster_hc_all_order_rename) <- paste("SO_level", 
        2:cluster_hc_dy_num, sep = "_")
    cluster_hc_all_order_rename_c <- apply(cluster_hc_all_order_rename, 
        2, function(x) {
            return(as.character(x))
        })
    row.names(cluster_hc_all_order_rename_c) <- row.names(cluster_hc_all_order_rename)
    cluster_hc_all_order_rename_c <- as.data.frame(cluster_hc_all_order_rename_c)
    for (i in 1:ncol(cluster_hc_all_order_rename_c)) {
        cluster_hc_all_order_rename_c[, i] <- paste("Term", i + 
            1, cluster_hc_all_order_rename_c[, i], sep = "_")
    }
    return(list("hclust_result_df" = cluster_hc_all_order_rename_c, 
        "hclust_result_model" = result_hc, "hclust_result_order" = cluster_hc_all_order))
}

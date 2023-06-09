DomainPlot <-
function (domain_hclust, distribution_distance, sample_information_coordinate=NULL, 
    k = ncol(domain_hclust[[3]]), size = 1, shape = 19, aspect_ratio = 1) 
{
    require(ggplot2)
    if(k>ncol(domain_hclust[[3]])){
        print(paste("The lagerest domain number is : ",ncol(domain_hclust[[3]]),"k cannot be larger!",sep=""))
              break
    }
    cluster_hc_all_order <- domain_hclust[[3]]
    cluster_hc = domain_hclust[[1]][, k]
    result_hc <- domain_hclust[[2]]
    names(cluster_hc) <- row.names(domain_hclust[[1]])
    a <- DrawCluster(t(distribution_distance), method = "umap", label = cluster_hc[row.names(distribution_distance)])
    colour_use <- ggplot_build(a$p)$data[[1]]$colour
    colour_use_unique <- unique(colour_use)
    names(colour_use_unique) <- 1:length(colour_use_unique)
    colour_use_unique_ok <- colour_use_unique[unique(cluster_hc_all_order[, 
        k])]
    ## A2Rplot draws hclust 
    source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    b <- A2Rplot(result_hc, k = k, boxes = FALSE, col.up = "gray50", 
        col.down = colour_use_unique_ok, show.labels = F, legend = T, 
        only.tree = T, type = c("triangle"))
    sample_information_coordinate$cluster_hc <- domain_hclust[[1]][row.names(sample_information_coordinate), 
        k]
    c <- ggplot(sample_information_coordinate, aes(x = X, y = Y, 
        colour = factor(cluster_hc))) + geom_point(size = size, 
        shape = shape) + theme(axis.text.x = element_text(vjust = 0.5, 
        hjust = 0.5, angle = 45), aspect.ratio = aspect_ratio)
    print(c)
}

DrawCluster <-
function (data, label = NULL, point_size = 1, method = c("tsne", 
    "umap"), draw_cluster_text = TRUE, calculated = TRUE, pca = TRUE, 
    perplexity = 100, plot = TRUE, legend_plot = TRUE, seed = 1) 
{
    require(ggplot2)
    require(dplyr)
    set.seed((seed))
    method = method[1]
    if (calculated) {
        if (method == "tsne") {
            require(Rtsne)
            tsneresult2 <- Rtsne(t(data), perplexity = perplexity, 
                pca = pca)
            X <- as.data.frame(tsneresult2$Y)
        }
        else if (method == "umap") {
            require(umap)
            umapresult1 <- umap(t(data))
            X <- as.data.frame(umapresult1$layout)
        }
        else {
            print("method must be tsne or umap.")
            break
        }
    }
    else {
        X = data
    }
    if (length(label) == 0) {
        label <- array(1, dim(X)[1])
        labelname = c(1)
    }
    labelname <- names(table(label))
    p <- ggplot(X, aes(x = X[, 1], y = X[, 2]))
    cell_group = factor(label)
    if (method == "tsne") {
        p <- p + geom_point(aes(color = cell_group), size = point_size) + 
            xlab("tSNE1") + ylab("tSNE2")
    }
    else {
        p <- p + geom_point(aes(color = cell_group), size = point_size) + 
            xlab("umap1") + ylab("umap2")
    }
    if (draw_cluster_text) {
        Label_cal <- X
        Label_cal$cluster <- label
        cluster_x_y <- Label_cal %>% group_by(cluster) %>% summarise(x_median = median(V1), 
            y_median = median(V2))
        p <- p + annotate("text", x = cluster_x_y$x_median, y = cluster_x_y$y_median, 
            label = cluster_x_y$cluster)
    }
    if (plot) {
        if (legend_plot) {
            mytheme <- theme_bw() + theme(plot.title = element_text(size = rel(1.5), 
                hjust = 0.5), axis.title = element_text(size = rel(1)), 
                axis.text = element_text(size = rel(1)), panel.grid.major = element_line(color = "white"), 
                panel.grid.minor = element_line(color = "white"), 
                legend.text = element_text(size = 10), legend.title = element_text(size = 15))
            p <- p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
            print(p)
        }
        else {
            mytheme <- theme_bw() + theme(plot.title = element_text(size = rel(1.5), 
                hjust = 0.5), axis.title = element_text(size = rel(1)), 
                axis.text = element_text(size = rel(1)), panel.grid.major = element_line(color = "white"), 
                panel.grid.minor = element_line(color = "white"), 
                legend.text = element_text(size = 10), legend.title = element_text(size = 15), 
                legend.position = "none")
            p <- p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
            print(p)
        }
    }
    return(list(p = p, x = X, cell_group = cell_group))
}

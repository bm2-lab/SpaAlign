DomainCharacterization <-
function (domain_hclust, cell_type_distribution,  
    k = ncol(domain_hclust[[1]]), single_cell = TRUE, sample_information_cellType=NULL,plot = TRUE) 
{
    cluster_hc = domain_hclust[[1]][, k]
    names(cluster_hc) <- row.names(domain_hclust[[1]])
    if(!single_cell){
        selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[1]])
        selected_propor <- cell_type_distribution[selected_samples, ]
        propor_mean <- apply(selected_propor, 2, mean)
        df <- data.frame(propor_mean)
        for (i in 2:length(table(cluster_hc))) {
            selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[i]])
            selected_propor <- cell_type_distribution[selected_samples, ]
            propor_mean <- apply(selected_propor, 2, mean)
            df[, i] <- propor_mean
        }
        colnames(df) <- names(table(cluster_hc))
        df_plot <- melt(df)
        df_plot$cellType <- rep(row.names(df), ncol(df))
        colnames(df_plot) <- c("Domain", "Proportion", "Cell_type")
        plot_c <- ggplot(df_plot, aes(Domain, Proportion, fill = Cell_type)) + 
            geom_bar(stat = "identity",color="white",size=0.2) + guides(fill = guide_legend(reverse = F)) + 
            scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(vjust = 0.5, 
            hjust = 0.5, angle = 45), aspect.ratio = 1)
        if(plot){
            plot(plot_c)
        }
        return(list(domain_average_distribution = df_plot))
    }else{
        sample_information_test<-sample_information_cellType
        cell_type_distribution_df <- as.data.frame(cell_type_distribution)
        cell_type_distribution_df$domain <- cluster_hc[row.names(cell_type_distribution_df)]
        selected_samples <- names(sample_information_test[sample_information_test == names(table(sample_information_test))[1]])
        selected_propor <- table(cluster_hc[selected_samples])/length(cluster_hc[selected_samples])
        df_cellType <- as.data.frame(matrix(0, length(table(sample_information_test)), length(table(cluster_hc))))
        row.names(df_cellType) <- names(table(sample_information_test))
        colnames(df_cellType) <- names(table(cluster_hc))
        df_cellType[1, names(selected_propor)] <- selected_propor
        for (i in 2:length(table(sample_information_test))) {
            selected_samples <- names(sample_information_test[sample_information_test == 
            names(table(sample_information_test))[i]])
            selected_propor <- table(cluster_hc[selected_samples])/length(cluster_hc[selected_samples])
            df_cellType[i, names(selected_propor)] <- selected_propor
        }
        occurence_cellType <- Correlation(t(df_cellType))
        
        selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[1]])
        selected_propor <- cell_type_distribution[selected_samples, ]
        propor_mean <- apply(selected_propor, 2, mean)
        df <- data.frame(propor_mean)
        for (i in 2:length(table(cluster_hc))) {
            selected_samples <- names(cluster_hc[cluster_hc == names(table(cluster_hc))[i]])
            selected_propor <- cell_type_distribution[selected_samples, ]
            propor_mean <- apply(selected_propor, 2, mean)
            df[, i] <- propor_mean
        }
        colnames(df) <- names(table(cluster_hc))
        df_plot <- melt(df)
        df_plot$cellType <- rep(row.names(df), ncol(df))
        colnames(df_plot) <- c("Domain", "Proportion", "Cell_type")
        
        
        df_cellType_plot <- melt(t(df_cellType))
        colnames(df_cellType_plot) <- c("Domain", "Cell_type", "Proportion")

        if (plot) {
            library(pheatmap)
            library(ggplot2)
            plot_b <- ggplot(df_cellType_plot, aes(Cell_type, Proportion, 
                fill = Domain)) + geom_bar(stat = "identity",color="white",size=0.2) + guides(fill = guide_legend(reverse = F)) + 
                scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(vjust = 0.5, 
                hjust = 0.5, angle = 45), aspect.ratio = 1)
            plot_c <- ggplot(df_plot, aes(Domain, Proportion, fill = Cell_type)) + 
                geom_bar(stat = "identity",color="white",size=0.2) + guides(fill = guide_legend(reverse = F)) + 
                scale_y_continuous(expand = c(0, 0)) + theme(axis.text.x = element_text(vjust = 0.5, 
                hjust = 0.5, angle = 45), aspect.ratio = 1)
            pheatmap(occurence_cellType)
            plot(plot_b)
            plot(plot_c)
        }
        return(list(occurence_rate = occurence_cellType, cell_type_average_distribution = df_cellType_plot, domain_average_distribution = df_plot))
    }
}

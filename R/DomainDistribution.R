DomainDistribution <-
function (hclust_result_df, datasets_label_combine, domain_num = ncol(hclust_result_df)) 
{
    domain_result_each_list <- list()
    SO_tree <- hclust_result_df[, domain_num]
    names(SO_tree) <- row.names(hclust_result_df)
    domain_names <- names(table(SO_tree))
    datasets_cellNum <- table(datasets_label_combine)
    datasets_num <- length(datasets_cellNum)
    datasets_name <- names(datasets_cellNum)
    proportion_dataset_domain <- matrix(0, domain_num, datasets_num)
    Term_table <- table(SO_tree)
    row.names(proportion_dataset_domain) <- names(Term_table)
    colnames(proportion_dataset_domain) <- names(datasets_cellNum)
    for (j in 1:length(Term_table)) {
        Datasets_term_cellNum <- table(datasets_label_combine[names(SO_tree[SO_tree == 
            names(Term_table)[j]])])
        proportion_dataset_domain[j, ] <- Datasets_term_cellNum/datasets_cellNum
    }
    for (i in 1:datasets_num) {
        domain_result_each_list[[i]] <- hclust_result_df[names(datasets_label_combine[datasets_label_combine == 
            datasets_name[i]]), ]
    }
    names(domain_result_each_list) <- datasets_name
    return(list("domain_distribution"=proportion_dataset_domain, "domain_result_each_list"=domain_result_each_list))
}

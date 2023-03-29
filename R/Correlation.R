Correlation <-
function (matrix, matrix2 = NULL, method = c("pearson", "spearman", 
    "cosine", "euclidean"), parallel_run=FALSE,cpu_num = ifelse(parallel_run==TRUE,8,FALSE)) 
{
    cosdist <- function(x1, x2) {
        n1 <- sqrt(sum(x1^2))
        n2 <- sqrt(sum(x2^2))
        d <- as.numeric(x1 %*% x2)/n1/n2
        return(d)
    }
    euclidean <- function(x1, x2) {
        return(sqrt(t(x1 - x2) %*% (x1 - x2)))
    }
    method <- method[1]
    if (!(method %in% c("pearson", "spearman", "cosine", "euclidean"))) {
        print("method must be 'pearson','spearman','cosine' or 'euclidean'.")
        break
    }
    simi_cor <- function(vec, matrix, method) {
        if (method == "cosine") {
            cor_matrix <- apply(matrix, 2, cosdist, x2 = vec)
        }
        else if (method == "euclidean") {
            cor_matrix <- apply(matrix, 2, euclidean, x2 = vec)
        }
        else {
            cor_matrix <- apply(matrix, 2, cor, y = vec, method = method)
        }
        return(cor_matrix)
    }
    sample_num <- ncol(matrix)
    cor_matrix <- matrix(rep(1, sample_num^2), sample_num)
    colnames(cor_matrix) <- colnames(matrix)
    row.names(cor_matrix) <- colnames(matrix)
    if(parallel_run){
        require(parallel)
        cpu_num_set <- makeCluster(getOption("cluster.cores", cpu_num))
        if (is.null(matrix2)) {
            cor_matrix <- parApply(cpu_num_set, matrix, 2, simi_cor, matrix = matrix, method = method)
            stopCluster(cpu_num_set)
            return(cor_matrix)
        }else {
            cor_matrix <- parApply(cpu_num_set, matrix, 2, simi_cor, matrix = matrix2, method = method)
            stopCluster(cpu_num_set)
            return(cor_matrix)
        }   
    }else{
        if (is.null(matrix2)) {
            cor_matrix <- apply(matrix, 2, simi_cor, matrix = matrix, method = method)
            return(cor_matrix)
        }else {
            cor_matrix <- apply(matrix, 2, simi_cor, matrix = matrix2, method = method)
            return(cor_matrix)
        }  
    }
}

`build.mim` <- function(data_matrix)
{
    data_matrix <- as.matrix(data_matrix)
    mi_matrix <- .Call("build_mim", as.vector(data_matrix), as.integer(nrow(data_matrix)),
            as.integer(ncol(data_matrix)), "ensemble");
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data_matrix), ncol=ncol(data_matrix))
    rownames(mi_matrix) <- colnames(data_matrix)
    colnames(mi_matrix) <- colnames(data_matrix)
    return(mi_matrix)
}

`filter.mRMR` <- function(children_count_per_level, feature_information_matrix,
        target_feature_index)
{
    children_count_per_level <- as.vector(children_count_per_level)
    feature_information_matrix <- as.vector(feature_information_matrix)
    target_feature_index <- as.integer(target_feature_index) - 1
    
    paths <- .Call("mRMR_filter", children_count_per_level, feature_information_matrix,
            target_feature_index, "ensemble")
    paths <- t(matrix(paths[length(paths):1], nrow=length(children_count_per_level),
                    ncol=length(paths)/length(children_count_per_level)))
    return(paths + 1)
}
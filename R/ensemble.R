`build.mim` <- function(data_matrix, strata=NULL, weights=NULL, feature_type=NULL)
{
	if(is.null(statra))
		strata <- rep.int(0, nrow(data_matrix))
		
	if(is.null(weights))
		weights <- rep.int(1, nrow(data_matrix))
		
	if(is.null(feature_type))
		feature_type <- rep.int(0, ncol(data_matrix))	
	
    data_matrix <- as.matrix(data_matrix)
    mi_matrix <- .Call("build_mim", as.vector(data_matrix), as.vector(strata), as.vector(weights), as.vector(feature_type), as.integer(nrow(data_matrix)),
            as.integer(ncol(data_matrix)), "ensemble");
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data_matrix), ncol=ncol(data_matrix))
    rownames(mi_matrix) <- colnames(data_matrix)
    colnames(mi_matrix) <- colnames(data_matrix)
    return(mi_matrix)
}

`filter.mRMR` <- function(data_matrix=NULL, feature_information_matrix=NULL,
        children_count_per_level, target_feature_index)
{
    children_count_per_level <- as.vector(children_count_per_level)
    target_feature_index <- as.integer(target_feature_index) - 1
    paths <- NULL
    
    if (!is.null(data_matrix))
    {
        data_matrix <- as.matrix(data_matrix)
        paths <- .Call("mRMR_filter_with_data", as.vector(data_matrix), as.integer(nrow(data_matrix)),
                as.integer(ncol(data_matrix)), children_count_per_level, feature_information_matrix,
                target_feature_index, "ensemble")
    }
    if (!is.null(feature_information_matrix))
    {
        paths <- .Call("mRMR_filter_with_mim", children_count_per_level,
                as.vector(feature_information_matrix), target_feature_index, "ensemble")
    }

    paths <- t(matrix(paths[length(paths):1], nrow=length(children_count_per_level),
                    ncol=length(paths)/length(children_count_per_level)))
    return(paths + 1)
}
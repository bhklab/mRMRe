`build.mim` <- function(data, strata=NULL, weights=NULL, feature_types=NULL)
{
    data <- as.matrix(data)
    
	if (is.null(strata))
		strata <- rep.int(0, nrow(data))
		
	if (is.null(weights))
		weights <- rep.int(1, nrow(data))
		
	if (is.null(feature_types))
		feature_types <- rep.int(0, ncol(data))	
	
    mi_matrix <- .Call("build_mim", as.vector(data), as.vector(strata), as.vector(weights),
            as.vector(feature_types), as.integer(nrow(data)), as.integer(ncol(data)),
            as.integer(length(unique(strata))), "ensemble");
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    rownames(mi_matrix) <- colnames(data)
    colnames(mi_matrix) <- colnames(data)
    
    return(mi_matrix)
}

`filter.mRMR` <- function(levels=NULL, data=NULL, strata=NULL, weights=NULL, feature_types=NULL,
        mim=NULL, target_feature_index=NULL)
{
    if (is.null(mim)) # filter_mRMR_with_data
    {
        data <- as.matrix(data)
        
        if (is.null(strata))
            strata <- rep.int(0, nrow(data))
        
        if (is.null(weights))
            weights <- rep.int(1, nrow(data))
        
        if (is.null(feature_types))
            feature_types <- rep.int(0, ncol(data))
        
        mRMR_tree <- .Call("filter_mRMR_with_data", as.vector(levels), data, as.vector(strata),
                as.vector(weights), as.vector(feature_types), nrow(data), ncol(data),
                as.integer(length(unique(strata))), as.integer(target_feature_index) - 1)
    }
    else if (is.null(data)) # filter_mRMR_with_mim
    {
        mim <- as.matrix(mim)
        
        mRMR_tree <- .Call("filter_mRMR_with_mim", as.vector(levels), mim, ncol(mim),
                as.integer(target_feature_index) - 1)   	
		
    }
	
	mRMR_tree$scores <- t(matrix(mRMR_tree$scores[length(mRMR_tree$scores):1], nrow=length(levels),
					ncol=length(mRMR_tree$scores)/length(levels)))
	mRMR_tree$paths <- t(matrix(mRMR_tree$paths[length(mRMR_tree$paths):1], nrow=length(levels),
					ncol=length(mRMR_tree$paths)/length(levels))) + 1

	return (mRMR_tree)
}

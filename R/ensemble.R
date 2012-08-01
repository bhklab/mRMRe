`build.mim` <- function(data, strata=NULL, weights=NULL, feature_types=NULL)
{
    data <- as.matrix(data)
    
	if (is.null(strata))
		strata <- rep.int(0, nrow(data))
		
	if (is.null(weights))
		weights <- rep.int(1, nrow(data))
		
	if (is.null(feature_types))
		feature_types <- rep.int(0, ncol(data))	
	
    mi_matrix <- .Call(build_mim, as.vector(data), as.vector(strata), as.vector(weights), as.vector(feature_types),
            as.integer(nrow(data)), as.integer(ncol(data)), as.integer(length(unique(strata))));
    mi_matrix <- matrix(mi_matrix, nrow=ncol(data), ncol=ncol(data))
    rownames(mi_matrix) <- colnames(data)
    colnames(mi_matrix) <- colnames(data)
    
    return(mi_matrix)
}

`filter.mRMR` <- function(levels, data=NULL, strata=NULL, weights=NULL, feature_types=NULL, mim=NULL,
        target_feature_index=NULL)
{
    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    if (is.null(mim))
    {
        data <- as.matrix(data)
        
        if (is.null(strata))
            strata <- rep.int(0, nrow(data))
        
        if (is.null(weights))
            weights <- rep.int(1, nrow(data))
        
        if (is.null(feature_types))
            feature_types <- rep.int(0, ncol(data))
        
        mRMR_tree <- .Call(filter_mRMR_with_data, as.vector(levels), data, as.vector(strata), as.vector(weights),
                as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
                as.integer(target_feature_index) - 1)
    }
    else if (is.null(data))
    {
        mim <- as.matrix(mim)
        
        mRMR_tree <- .Call(filter_mRMR_with_mim, as.vector(levels), mim, ncol(mim),
                as.integer(target_feature_index) - 1)
    }
	
    mRMR_tree$paths <- wrap(mRMR_tree$paths) + 1
	mRMR_tree$scores <- wrap(mRMR_tree$scores)
    
	return (mRMR_tree)
}
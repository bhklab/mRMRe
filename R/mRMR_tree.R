`mRMR_tree` <- function(paths, scores)
{
    object <- list(paths=paths, scores=scores)
    class(object) <- "mRMR_tree"
    return(object)
}

`filter.mRMR_tree` <- function(levels, data=NULL, strata=NULL, weights=NULL, feature_types=NULL, mim=NULL,
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
        
        tree <- .Call(C_build_mRMR_tree_from_data, as.vector(levels), data, as.vector(strata), as.vector(weights),
                as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
                as.integer(target_feature_index) - 1)
    }
    else if (is.null(data))
    {
        mim <- as.matrix(mim)
        
        tree <- .Call(C_build_mRMR_tree_from_mim, as.vector(levels), mim, ncol(mim),
                as.integer(target_feature_index) - 1)
    }
    
    return(mRMR_tree(paths=wrap(mRMR_tree$paths) + 1, scores=wrap(mRMR_tree$scores)))
}
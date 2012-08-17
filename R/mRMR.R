`mRMR.classic` <- function(
        data,
        target,
        feature_count,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=rep.int(0, ncol(data) + 1), # First one must be the target
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0
        )
{
    return(mRMRe::mRMR.ensemble(levels=rep.int(1, feature_count), data=data, target=target, strata=strata,
                    weights=weights, feature_types=feature_types, uses_ranks=uses_ranks, outX=outX,
                    bootstrap_count=bootstrap_count))
}

# feature_type = 0 for continous
# feature_type = 1 for discrete
# feature_type = 2 for survival
`mRMR.ensemble` <- function(
        data,
        target,
        levels,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=rep.int(0, ncol(data) + 1), # First one must be the target
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    if (is.data.frame(data))
        stop("data must be of type data frame")
    
    levels <- as.vector(levels)
    data <- cbind(target=target, data)

    expansion <- expand.data(data, feature_types)
    data <- expansion$data
    feature_types <- expansion$feature_types
    
    tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), as.vector(strata), as.vector(weights),
            as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
            0, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
	tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))
    rownames(tree$mim) <- colnames(data)
    colnames(tree$mim) <- colnames(data)

    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    object <- list(paths=wrap(tree$paths) + 1, scores=wrap(tree$scores), tree$mim)
    class(object) <- "mRMReObject"
    return(object)
}
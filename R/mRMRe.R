`mRMRclassic` <- function(
        data,
        target,
        feature_count,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=rep.int(0, ncol(data)),
        mim=NULL,
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0
        )
{
    return(mRMRe::mRMRensemble(levels=rep.int(1, feature_count), data=data, target=target, strata=strata,
                    weights=weights, feature_types=feature_types, mim=mim, uses_ranks=uses_ranks, outX=outX,
                    bootstrap_count=bootstrap_count))
}

`mRMRensemble` <- function(
        data=NULL,
        target=NULL,
        levels,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=rep.int(0, ncol(data)),
        mim=NULL,
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    levels <- as.vector(levels)
    data <- cbind(target, data)
    
    if (is.null(mim))
    {
        tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), as.vector(strata), as.vector(weights),
                as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
                0, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
		tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))
        rownames(tree$mim) <- colnames(data)
        colnames(tree$mim) <- colnames(data)
    }
    else if (is.null(data))
    {
        tree <- .Call(C_build_mRMR_tree_from_mim, levels, as.vector(mim), ncol(mim),
                target_feature_index)
        tree$mim <- mim
    }

    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    object <- list(paths=wrap(tree$paths) + 1, scores=wrap(tree$scores), tree$mim)
    class(object) <- "mRMRensemble"
    return(object)
}

`set_thread_count` <- function(thread_count)
{
    .Call(C_set_thread_count, as.integer(thread_count))
}


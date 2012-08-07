`mRMR_tree` <- function(
        paths,
        scores,
        target_feature_index,
        mim)
{
    object <- list(paths=paths, scores=scores, target_feature_index=target_feature_index, mim=mim)
    class(object) <- "mRMR_tree"
    return(object)
}

`filter.mRMR_tree` <- function(
        levels,
        data=NULL,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=rep.int(0, ncol(data)),
        mim=NULL,
        target_feature_index=NULL,
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    levels <- as.vector(levels)
    target_feature_index <- as.integer(target_feature_index) - 1
    
    if (is.null(mim))
    {
        tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), as.vector(strata), as.vector(weights),
                as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
                target_feature_index, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
    }
    else if (is.null(data))
    {
        tree <- .Call(C_build_mRMR_tree_from_mim, levels, as.vector(mim), ncol(mim),
                target_feature_index)
        tree$mim <- mim
    }

    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))
    return(mRMR_tree(paths=wrap(tree$paths) + 1, scores=wrap(tree$scores), target_feature_index + 1, tree$mim))
}
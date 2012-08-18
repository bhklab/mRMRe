`mRMR.classic` <- function(
        data,
        target,
        feature_count,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0
        )
{
    return(mRMRe::mRMR.ensemble(levels=rep.int(1, feature_count), data=data, target=target, strata=strata,
                    weights=weights, uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count))
}

`mRMR.ensemble` <- function(
        data,
        target_index,
        levels,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    if (is.data.frame(data))
        stop("data must be of type data frame") }
    
    feature_types <- unlist(sapply(data, class)
    if (is.list(feature_types))
        feature_types <- unlist(lapply(feature_types, paste, collapse="_"))

    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("feature types must be either numeric, ordered factor or Surv")

    # if (any(!sapply(data[ ,feature_types == "factor", drop=FALSE], is.ordered)))
    #     stop("categorical features must be ordered factors")

    expansion <- .expand.data(data=data)
    data <- expansion$data
    feature_types <- expansion$feature_types
    levels <- as.vector(levels)
        
    tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), as.vector(strata), as.vector(weights),
            as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
            as.integer(target_index) - 1, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
    tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))
    rownames(tree$mim) <- colnames(data)
    colnames(tree$mim) <- colnames(data)

    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    paths <- wrap(tree$paths)
    
    offsets <- rep(0, length(paths))
    lapply(which(feature_types == 3), function(index)
    {
        subset <- which(paths > index)
        offsets[subset] <<- offsets[subset] + 1
    })
    paths <- paths - offsets
    
    object <- list(target_index=target_index, paths=paths, scores=wrap(tree$scores), mim=tree$mim)
    class(object) <- "mRMReObject"
    return(object)
}

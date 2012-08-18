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
        target_index,
        levels,
        strata=rep.int(0, nrow(data)),
        weights=rep.int(1, nrow(data)),
        feature_types=c("numeric", "factor", "Surv"), # First one must be the target
        uses_ranks=TRUE,
        outX=TRUE,
        bootstrap_count=0)
{
    if (is.data.frame(data)) { stop("data must be of type data frame") }
    feature_types <- sapply(data, class)
    if(is.list(feature_types)) { feature_types <- unlist(lapply(feature_types, function(x) { paste(x, collapse="_") })) }
    if(any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv")))) { stop("feature types must be either numeric, ordered factor or Surv") }
    #if(any(!sapply(data[ ,feature_types == "factor", drop=FALSE], is.ordered))) { stop("categorical features must be ordered factors") }
    expansion <- .expand.data(data=data)
    data <- expansion$data
    feature_types <- expansion$feature_types

    levels <- as.vector(levels)
        
    tree <- .Call(C_build_mRMR_tree_from_data, levels, as.vector(data), as.vector(strata), as.vector(weights),
            as.vector(feature_types), nrow(data), ncol(data), as.integer(length(unique(strata))),
            0, as.integer(uses_ranks), as.integer(outX), as.integer(bootstrap_count))
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
    
    object <- list(paths=paths, scores=wrap(tree$scores), mim=tree$mim)
    class(object) <- "mRMReObject"
    return(object)
}
`mRMR.classic` <- function(data, priors, prior_weights, target_index, feature_count, strata, weights, uses_ranks,
        outX, bootstrap_count)
{
    return(mRMRe::mRMR.filter(levels=rep.int(1, feature_count), data=data, priors=priors,
                    prior_weights=prior_weights, target=target_index, strata=strata, weights=weights,
                    uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count))
}

`mRMR.ensemble` <- function(data, priors, prior_weights, target_index, feature_count, solution_count, strata, weights,
        uses_ranks, outX, bootstrap_count)
{
    return(mRMRe::mRMR.filter(levels=c(solution_count, rep.int(1, feature_count - 1)), data=data, priors=priors,
                    prior_weights=prior_weights, target=target_index, strata=strata, weights=weights,
                    uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count))
}

`mRMR.filter` <- function(data, priors, prior_weights, target_index, levels, strata, weights, uses_ranks=TRUE,
        outX=TRUE, bootstrap_count=0)
{
    if (!is.data.frame(data))
        stop("data must be of type data frame")
    
    feature_types <- sapply(data, function(x) paste(class(x), collapse="_"))
    
    if (any(!is.element(feature_types, c("numeric", "ordered_factor", "Surv"))))
        stop("data columns must be either of numeric, ordered factor or Surv type")

    if (missing(weights)) 
        weights <- rep.int(1, nrow(data))
    
    if (missing(strata)) 
        strata <- rep.int(0, nrow(data))
    else if (is.factor(strata))
        strata <- as.integer(strata) - 1
    else
        stop("strata must be provided as factors")
    
    if (missing(priors))
    {
        priors <- vector()
        prior_weights <- 0
    }
    else if (missing(prior_weights))
        stop("prior_weights must be provided with priors")
    else if (prior_weights > 1 || prior_weights < 0)
        stop("prior_weights must be [0, 1]")
	else if (max(priors) > 1 || min(priors) < 0)
		stop("prior must be [0,1]")
    
    levels <- as.vector(levels)
    expansion <- mRMRe:::.expand.input(feature_types=feature_types, data=data, priors=priors)
    data <- expansion$data
    feature_types <- expansion$feature_types
    feature_names <- expansion$feature_names
    priors <- expansion$priors
    
    tree <- .Call(mRMRe:::.C_build_mRMR_tree, levels, as.vector(data), as.vector(priors), as.numeric(prior_weights),
            as.vector(strata), as.vector(weights), as.vector(feature_types), nrow(data), ncol(data),
            as.integer(length(unique(strata))), as.integer(target_index) - 1, as.integer(uses_ranks),
            as.integer(outX), as.integer(bootstrap_count))
    tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))

    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    compression <- mRMRe:::.compress.output(feature_types=feature_types, feature_names=feature_names, mi_matrix=tree$mim, paths=wrap(tree$paths))
    mi_matrix <- compression$mi_matrix
    paths <- compression$paths + 1
    
    object <- list("target_index"=target_index, "paths"=paths, "scores"=wrap(tree$scores), "mi_matrix"=mi_matrix)
    class(object) <- "mRMReFilter"
    return(object)
}

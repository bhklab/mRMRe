`mRMR.classic` <- function(data, priors, prior_weights, target_index, feature_count, strata, weights, uses_ranks,
        outX, bootstrap_count)
{
    return(mRMRe::mRMR.filter(data=data, priors=priors, prior_weights=prior_weights, target=target_index,
                    levels=rep.int(1, feature_count), strata=strata, weights=weights, uses_ranks=uses_ranks,
                    outX=outX, bootstrap_count=bootstrap_count))
}

`mRMR.ensemble` <- function(data, priors, prior_weights, target_index, feature_count, solution_count, strata, weights,
        uses_ranks, outX, bootstrap_count)
{
    return(mRMRe::mRMR.filter(data=data, priors=priors, prior_weights=prior_weights, target=target_index,
                    levels=c(solution_count, rep.int(1, feature_count - 1)), strata=strata, weights=weights,
                    uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count))
}

`mRMR.filter` <- function(data, priors, prior_weights, target_index, levels, strata, weights, uses_ranks, outX,
        bootstrap_count)
{
    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    ## Expansion
    expansion <- mRMRe:::.expand.input(data=data, priors=priors, prior_weights=prior_weights, strata=strata,
            weights=weights)
    data <- expansion$data
    priors <- expansion$priors
    prior_weights <- expansion$prior_weights
    strata <- expansion$strata
    weights <- expansion$weights
    feature_types <- expansion$feature_types
    feature_names <- expansion$feature_names
    uses_ranks <- expansion$uses_ranks
    outX <- expansion$outX
    bootstrap_count <- expansion$bootstrap_count
    
    ## Call
    tree <- .Call(mRMRe:::.C_build_mRMR_tree, as.vector(levels), as.vector(data), as.vector(priors),
            as.numeric(prior_weights), as.vector(strata), as.vector(weights), as.vector(feature_types), nrow(data),
            ncol(data), as.integer(length(unique(strata))), as.integer(target_index) - 1, as.integer(uses_ranks),
            as.integer(outX), as.integer(bootstrap_count))
    tree$mim <- matrix(tree$mim, ncol=sqrt(length(tree$mim)), nrow=sqrt(length(tree$mim)))
    tree$paths <- wrap(tree$paths)
    tree$scores <- wrap(tree$scores)
    
    ## Compression
    compression <- mRMRe:::.compress.output(feature_types=feature_types, feature_names=feature_names,
            mi_matrix=tree$mim, paths=tree$paths)
    mi_matrix <- compression$mi_matrix
    paths <- compression$paths + 1
    
    ## Return
    object <- list("target_index"=target_index, "paths"=paths, "scores"=tree$scores, "mi_matrix"=mi_matrix)
    class(object) <- "mRMReFilter"
    return(object)
}

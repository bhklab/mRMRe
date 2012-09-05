`mRMR.classic` <- function(data, priors, prior_weight, target_index, feature_count, strata, weights, uses_ranks,
        outX, bootstrap_count, .is_expanded, .feature_types, .feature_names)
{
    return(mRMRe::mRMR.ensemble(data=data, priors=priors, prior_weight=prior_weight, target=target_index,
                    feature_count=feature_count, solution_count=1, strata=strata, weights=weights,
                    uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count, .is_expanded=.is_expanded,
                    .feature_types=.feature_types, .feature_names=.feature_names))
}

`mRMR.ensemble` <- function(data, priors, prior_weight, target_index, feature_count, solution_count, strata, weights,
        uses_ranks, outX, bootstrap_count, .is_expanded, .feature_types, .feature_names)
{
    return(mRMRe::mRMR.filter(data=data, priors=priors, prior_weight=prior_weight, target=target_index,
                    levels=c(solution_count, rep.int(1, feature_count - 1)), strata=strata, weights=weights,
                    uses_ranks=uses_ranks, outX=outX, bootstrap_count=bootstrap_count, .is_expanded=.is_expanded,
                    .feature_types=.feature_types, .feature_names=.feature_names))
}

`mRMR.filter` <- function(data, priors, prior_weight, target_index, levels, strata, weights, uses_ranks, outX,
        bootstrap_count, .is_expanded, .feature_types, .feature_names)
{
    feature_names <- colnames(data)
    if (missing(.is_expanded) || !.is_expanded)
    {
        expansion <- mRMRe:::.expand.input(data=data, priors=priors, prior_weight=prior_weight, strata=strata,
                weights=weights, target_indices=target_index, uses_ranks=uses_ranks, outX=outX,
                bootstrap_count=bootstrap_count)
        data <- expansion$data
        priors <- expansion$priors
        prior_weight <- expansion$prior_weight
        strata <- expansion$strata
        weights <- expansion$weights
        feature_types <- expansion$feature_types
        uses_ranks <- expansion$uses_ranks
        outX <- expansion$outX
        bootstrap_count <- expansion$bootstrap_count
        target_index <- expansion$target_indices
    }
    else
    {
        if (!missing(.feature_types))
            feature_types <- .feature_types
        else
            stop(".feature_types must be provided when .is_expanded is set to TRUE")
        
        if (!missing(.feature_names))
            feature_names <- .feature_names
        else
            stop(".feature_names must be provided when .is_expanded is set to TRUE")
    }
    
    wrap <- function(i) t(matrix(i[length(i):1], nrow=length(levels), ncol=length(i)/length(levels)))
    
    tree <- .Call(mRMRe:::.C_export_filter, as.vector(levels), as.vector(data), as.vector(priors),
            as.numeric(prior_weight), as.vector(strata), as.vector(weights), as.vector(feature_types), nrow(data),
            ncol(data), as.integer(length(unique(strata))), as.integer(target_index) - 1, as.integer(uses_ranks),
            as.integer(outX), as.integer(bootstrap_count))
    tree$mi_matrix <- matrix(tree$mi_matrix, ncol=sqrt(length(tree$mi_matrix)), nrow=sqrt(length(tree$mi_matrix)))
    tree$paths <- wrap(tree$paths)

    compression <- mRMRe:::.compress.output(feature_types=feature_types, feature_names=feature_names,
            mi_matrix=tree$mi_matrix, paths=tree$paths, target_indices=target_index)
    mi_matrix <- compression$mi_matrix
    paths <- compression$paths + 1
    target_index <- compression$target_indices
    
    object <- list("target_index"=target_index, "paths"=paths, "mi_matrix"=mi_matrix)
    class(object) <- "mRMReFilter"
    return(object)
}

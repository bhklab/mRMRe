`mRMR.network` <- function(data, priors, prior_weights, feature_indices, feature_count, solution_count, strata,
        weights, uses_ranks, outX, bootstrap_count, layers)
{
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
    
    topologies <- list()
    length(topologies) <- ncol(data)
    
    lapply(seq(layers), function(layer)
    {
        feature_indices <<- intersect(unlist(lapply(feature_indices, function(target_index)
        {
            filter <- mRMRe::mRMR.ensemble(data=data, priors=priors, prior_weights=prior_weights,
                    target_index=target_index, feature_count=feature_count, solution_count=solution_count,
                    strata=strata, weights=weights, uses_ranks=uses_ranks, outX=outX,
                    bootstrap_count=bootstrap_count)
            topologies[[target_index]] <- filter$paths # Add scoring scheme here
            return(as.vector(filter$paths))
        })), which(!is.null(topologies)))
        return(NULL)
    })
    
    object <- list("topologies"=topologies)
    class(object) <- "mRMReNetwork"
    return(object)
}
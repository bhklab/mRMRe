`mRMR.network` <- function(data, priors, prior_weight, target_indices, feature_count, solution_count, strata,
        weights, uses_ranks, outX, bootstrap_count, layers)
{
    if (missing(layers))
        layers <- 1L

    topologies <- list()
    length(topologies) <- ncol(data)
    
    old_feature_types <- sapply(data, function(x) paste(class(x), collapse="_"))
    
    
    feature_names <- colnames(data)
    expansion <- mRMRe:::.expand.input(data=data, priors=priors, prior_weight=prior_weight, strata=strata,
            weights=weights, target_indices=target_indices)
    data <- expansion$data
    priors <- expansion$priors
    prior_weight <- expansion$prior_weight
    strata <- expansion$strata
    weights <- expansion$weights
    feature_types <- expansion$feature_types
    uses_ranks <- expansion$uses_ranks
    outX <- expansion$outX
    bootstrap_count <- expansion$bootstrap_count
    target_indices <- expansion$target_indices

    lapply(seq(layers), function(layer)
    {
        target_indices <<- unlist(lapply(target_indices, function(target_index)
        {
            filter <- mRMRe::mRMR.ensemble(data=data, priors=priors, prior_weight=prior_weight,
                    target_index=target_index, feature_count=feature_count, solution_count=solution_count,
                    strata=strata, weights=weights, uses_ranks=uses_ranks, outX=outX,
                    bootstrap_count=bootstrap_count, .is_expanded=TRUE, .feature_types=feature_types,
                    .feature_names=feature_names)
            topologies[[filter$target_index]] <<- filter$paths # Add scoring scheme here
            return(as.vector(filter$paths))
        }))
        target_indices <<- intersect(target_indices, which(sapply(topologies, is.null)))
        target_indices <<- mRMRe:::.expand.feature.indices(old_feature_types, target_indices)
        return(NULL)
    })
    
    object <- list("topologies"=topologies)
    class(object) <- "mRMReNetwork"
    return(object)
}